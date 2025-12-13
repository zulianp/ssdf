#include <mpi.h>
#include <algorithm>
#include <limits>
#include <numeric>
#include <vector>
#include "ssdf.hpp"

#define PSDF_TIMER(comm, name) ssdf::BarrierTimer t_##name(comm, #name);
// #define PSDF_TIMER(...)

namespace ssdf {

    struct BarrierTimer {
        MPI_Comm comm;
        std::string name;
        double start_time;
        double end_time;

        inline BarrierTimer(MPI_Comm comm, const std::string &name) : comm(comm), name(name) {
            MPI_Barrier(comm);
            start_time = time_ms();
        }

        inline ~BarrierTimer() {
            end_time = time_ms();

            int rank, size;
            MPI_Comm_rank(comm, &rank);
            MPI_Comm_size(comm, &size);

            double times[2] = {start_time, end_time};

            if (!rank) {
                std::vector<double> all(size * 2);
                MPI_Gather(times, 2, MPI_DOUBLE, all.data(), 2, MPI_DOUBLE, 0, comm);

                double begin = all[0];
                double end = all[1];
                double avg = 0;

                for (int r = 0; r < size; r++) {
                    const int i = r * 2;
                    const double s = all[i];
                    const double e = all[i + 1];
                    avg += e - s;

                    begin = std::min(begin, s);
                    end = std::max(end, e);
                }
                avg /= size;

                printf("%s: overall %g [ms], avg %g [ms]\n", name.c_str(), end - begin, avg);
            } else {
                MPI_Gather(times, 2, MPI_DOUBLE, nullptr, 0, MPI_DOUBLE, 0, comm);
            }
        }
    };

    template <typename T>
    MPI_Datatype mpi_type();

    template <>
    MPI_Datatype mpi_type<float>() {
        return MPI_FLOAT;
    }
    template <>
    MPI_Datatype mpi_type<double>() {
        return MPI_DOUBLE;
    }
    template <>
    MPI_Datatype mpi_type<int>() {
        return MPI_INT;
    }
    template <>
    MPI_Datatype mpi_type<long long>() {
        return MPI_LONG_LONG;
    }

    /**
     * @brief Distributed unsigned point-to-surface distances with MPI.
     *
     * @tparam G Geometry type (float/double) for coordinates.
     * @tparam T Output distance type.
     * @tparam I Index type for triangle vertices.
     * @param comm MPI communicator.
     * @param lnpoints Local number of query points on this rank.
     * @param x,y,z Local point coordinates (SoA), size lnpoints.
     * @param lnselements Local number of surface triangles.
     * @param s0,s1,s2 Local triangle indices (size lnselements), referencing local surface vertices.
     * @param lnspoints Local number of surface vertices.
     * @param sx,sy,sz Local surface vertex coordinates (SoA), size lnspoints.
     * @param out Output array of size lnpoints. On input, values are treated as current best
     *            distances (initialize to large values if unused).
     * @return int 0 on success.
     *
     * Notes:
     * - Distances are unsigned (closest-point Euclidean distance).
     * - Data distribution is contiguous partitions of points and surface across ranks.
     * - OpenMP may parallelize local computation if enabled at build time.
     */
    template <typename G, typename T, typename I>
    int pedf(MPI_Comm comm,
             const ptrdiff_t lnpoints,
             const G *const SSDF_RESTRICT x,
             const G *const SSDF_RESTRICT y,
             const G *const SSDF_RESTRICT z,
             const ptrdiff_t lnselements,  // local number of surface elements
             const I *const SSDF_RESTRICT s0,
             const I *const SSDF_RESTRICT s1,
             const I *const SSDF_RESTRICT s2,
             const ptrdiff_t lnspoints,  // local number of surcace points
             const G *const SSDF_RESTRICT sx,
             const G *const SSDF_RESTRICT sy,
             const G *const SSDF_RESTRICT sz,
             T *const SSDF_RESTRICT out) {
        PSDF_TIMER(comm, total);
        int comm_size = 0, comm_rank = 0;
        MPI_Comm_size(comm, &comm_size);
        MPI_Comm_rank(comm, &comm_rank);

        auto mpi_type_G = mpi_type<G>();
        auto mpi_type_I = mpi_type<I>();

        int SSDF_USE_ALLGATHER = 1;
        SSDF_READ_ENV(SSDF_USE_ALLGATHER, atoi);

        // Ensure output starts with large values (DONE Outside)
        // for (ptrdiff_t i = 0; i < lnpoints; ++i) {
        //     out[i] = std::numeric_limits<T>::max();
        // }

        if (SSDF_USE_ALLGATHER) {
            // Gather counts
            std::vector<long long> all_spoints(comm_size), all_selems(comm_size);
            long long lnspoints_ll = static_cast<long long>(lnspoints);
            long long lnselements_ll = static_cast<long long>(lnselements);
            MPI_Allgather(&lnspoints_ll, 1, MPI_LONG_LONG, all_spoints.data(), 1, MPI_LONG_LONG, comm);
            MPI_Allgather(&lnselements_ll, 1, MPI_LONG_LONG, all_selems.data(), 1, MPI_LONG_LONG, comm);

            // Compute offsets
            std::vector<long long> off_spoints(comm_size + 1, 0);
            std::vector<long long> off_selems(comm_size + 1, 0);
            for (int r = 0; r < comm_size; ++r) {
                off_spoints[r + 1] = off_spoints[r] + all_spoints[r];
                off_selems[r + 1] = off_selems[r] + all_selems[r];
            }
            const ptrdiff_t nspoints_total = static_cast<ptrdiff_t>(off_spoints.back());
            const ptrdiff_t nselements_total = static_cast<ptrdiff_t>(off_selems.back());

            // Prepare local, index-adjusted copies
            std::vector<I> s0_adj(lnselements), s1_adj(lnselements), s2_adj(lnselements);
            const I point_offset = static_cast<I>(off_spoints[comm_rank]);
            for (ptrdiff_t i = 0; i < lnselements; ++i) {
                s0_adj[i] = s0[i] + point_offset;
                s1_adj[i] = s1[i] + point_offset;
                s2_adj[i] = s2[i] + point_offset;
            }

            // Allocate global buffers
            std::vector<G> gx(nspoints_total), gy(nspoints_total), gz(nspoints_total);
            std::vector<I> gs0(nselements_total), gs1(nselements_total), gs2(nselements_total);

            // Build counts/displs
            std::vector<int> sc_sp(comm_size), sc_se(comm_size), dsp_sp(comm_size), dsp_se(comm_size);
            for (int r = 0; r < comm_size; ++r) {
                sc_sp[r] = static_cast<int>(all_spoints[r]);
                sc_se[r] = static_cast<int>(all_selems[r]);
                dsp_sp[r] = static_cast<int>(off_spoints[r]);
                dsp_se[r] = static_cast<int>(off_selems[r]);
            }

            MPI_Allgatherv(sx, sc_sp[comm_rank], mpi_type_G, gx.data(), sc_sp.data(), dsp_sp.data(), mpi_type_G, comm);
            MPI_Allgatherv(sy, sc_sp[comm_rank], mpi_type_G, gy.data(), sc_sp.data(), dsp_sp.data(), mpi_type_G, comm);
            MPI_Allgatherv(sz, sc_sp[comm_rank], mpi_type_G, gz.data(), sc_sp.data(), dsp_sp.data(), mpi_type_G, comm);

            MPI_Allgatherv(
                s0_adj.data(), sc_se[comm_rank], mpi_type_I, gs0.data(), sc_se.data(), dsp_se.data(), mpi_type_I, comm);
            MPI_Allgatherv(
                s1_adj.data(), sc_se[comm_rank], mpi_type_I, gs1.data(), sc_se.data(), dsp_se.data(), mpi_type_I, comm);
            MPI_Allgatherv(
                s2_adj.data(), sc_se[comm_rank], mpi_type_I, gs2.data(), sc_se.data(), dsp_se.data(), mpi_type_I, comm);

            // Compute EDF for local points against global surface
            ssdf::edf<G, T, I>(lnpoints,
                               x,
                               y,
                               z,
                               nselements_total,
                               gs0.data(),
                               gs1.data(),
                               gs2.data(),
                               nspoints_total,
                               gx.data(),
                               gy.data(),
                               gz.data(),
                               out);
        } else {
            // Reduce max sizes to size buffers
            long long max_spoints = static_cast<long long>(lnspoints);
            long long max_selems = static_cast<long long>(lnselements);
            MPI_Allreduce(MPI_IN_PLACE, &max_spoints, 1, MPI_LONG_LONG, MPI_MAX, comm);
            MPI_Allreduce(MPI_IN_PLACE, &max_selems, 1, MPI_LONG_LONG, MPI_MAX, comm);

            std::vector<G> sx_buf(max_spoints), sy_buf(max_spoints), sz_buf(max_spoints);
            std::vector<I> s0_buf(max_selems), s1_buf(max_selems), s2_buf(max_selems);

            for (int root = 0; root < comm_size; ++root) {
                PSDF_TIMER(comm, round);

                long long nsp_root = 0, nse_root = 0;
                if (comm_rank == root) {
                    nsp_root = lnspoints;
                    nse_root = lnselements;
                }
                MPI_Bcast(&nsp_root, 1, MPI_LONG_LONG, root, comm);
                MPI_Bcast(&nse_root, 1, MPI_LONG_LONG, root, comm);

                if (nsp_root == 0 || nse_root == 0) continue;

                if (comm_rank == root) {
                    std::copy(sx, sx + nsp_root, sx_buf.begin());
                    std::copy(sy, sy + nsp_root, sy_buf.begin());
                    std::copy(sz, sz + nsp_root, sz_buf.begin());
                    std::copy(s0, s0 + nse_root, s0_buf.begin());
                    std::copy(s1, s1 + nse_root, s1_buf.begin());
                    std::copy(s2, s2 + nse_root, s2_buf.begin());
                }

                MPI_Bcast(sx_buf.data(), static_cast<int>(nsp_root), mpi_type_G, root, comm);
                MPI_Bcast(sy_buf.data(), static_cast<int>(nsp_root), mpi_type_G, root, comm);
                MPI_Bcast(sz_buf.data(), static_cast<int>(nsp_root), mpi_type_G, root, comm);
                MPI_Bcast(s0_buf.data(), static_cast<int>(nse_root), mpi_type_I, root, comm);
                MPI_Bcast(s1_buf.data(), static_cast<int>(nse_root), mpi_type_I, root, comm);
                MPI_Bcast(s2_buf.data(), static_cast<int>(nse_root), mpi_type_I, root, comm);

                // Compute EDF against this chunk
                ssdf::edf<G, T, I>(lnpoints,
                                   x,
                                   y,
                                   z,
                                   static_cast<ptrdiff_t>(nse_root),
                                   s0_buf.data(),
                                   s1_buf.data(),
                                   s2_buf.data(),
                                   static_cast<ptrdiff_t>(nsp_root),
                                   sx_buf.data(),
                                   sy_buf.data(),
                                   sz_buf.data(),
                                   out);
            }
        }

        return 0;
    }
}  // namespace ssdf