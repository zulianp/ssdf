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
            start_time = MPI_Wtime();
        }

        inline ~BarrierTimer() {
            end_time = MPI_Wtime();

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

                printf("%s: overall %g [s], avg %g [s]\n", name.c_str(), end - begin, avg);
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
            ssdf::edf_select<G, T, I>(lnpoints,
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
                ssdf::edf_select<G, T, I>(lnpoints,
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


    template <typename G, typename T, typename I>
    int pedf_update(MPI_Comm comm,
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
            PSDF_TIMER(comm, pedf_update);
            int comm_size = 0, comm_rank = 0;
            MPI_Comm_size(comm, &comm_size);
            MPI_Comm_rank(comm, &comm_rank);

            auto mpi_type_G = mpi_type<G>();
            auto mpi_type_I = mpi_type<I>();

            // 1) Single bounding box per rank around local points inflated by current distances
            G pminx = std::numeric_limits<G>::max(), pminy = std::numeric_limits<G>::max(),
              pminz = std::numeric_limits<G>::max();
            G pmaxx = std::numeric_limits<G>::lowest(), pmaxy = std::numeric_limits<G>::lowest(),
              pmaxz = std::numeric_limits<G>::lowest();
#pragma omp parallel
            {
                G tminx = std::numeric_limits<G>::max(), tminy = std::numeric_limits<G>::max(),
                  tminz = std::numeric_limits<G>::max();
                G tmaxx = std::numeric_limits<G>::lowest(), tmaxy = std::numeric_limits<G>::lowest(),
                  tmaxz = std::numeric_limits<G>::lowest();
#pragma omp for nowait
                for (ptrdiff_t i = 0; i < lnpoints; ++i) {
                    const G r = static_cast<G>(out[i]);
                    const G lx = x[i], ly = y[i], lz = z[i];
                    tminx = std::min(tminx, lx - r);
                    tmaxx = std::max(tmaxx, lx + r);
                    tminy = std::min(tminy, ly - r);
                    tmaxy = std::max(tmaxy, ly + r);
                    tminz = std::min(tminz, lz - r);
                    tmaxz = std::max(tmaxz, lz + r);
                }
#pragma omp critical
                {
                    pminx = std::min(pminx, tminx);
                    pmaxx = std::max(pmaxx, tmaxx);
                    pminy = std::min(pminy, tminy);
                    pmaxy = std::max(pmaxy, tmaxy);
                    pminz = std::min(pminz, tminz);
                    pmaxz = std::max(pmaxz, tmaxz);
                }
            }

            // 2) Allgather one bounding box per rank (fixed-size)
            std::vector<G> all_pminx(comm_size), all_pmaxx(comm_size), all_pminy(comm_size), all_pmaxy(comm_size),
                all_pminz(comm_size), all_pmaxz(comm_size);
            MPI_Allgather(&pminx, 1, mpi_type_G, all_pminx.data(), 1, mpi_type_G, comm);
            MPI_Allgather(&pmaxx, 1, mpi_type_G, all_pmaxx.data(), 1, mpi_type_G, comm);
            MPI_Allgather(&pminy, 1, mpi_type_G, all_pminy.data(), 1, mpi_type_G, comm);
            MPI_Allgather(&pmaxy, 1, mpi_type_G, all_pmaxy.data(), 1, mpi_type_G, comm);
            MPI_Allgather(&pminz, 1, mpi_type_G, all_pminz.data(), 1, mpi_type_G, comm);
            MPI_Allgather(&pmaxz, 1, mpi_type_G, all_pmaxz.data(), 1, mpi_type_G, comm);

            auto intersects = [](G minx1, G maxx1, G miny1, G maxy1, G minz1, G maxz1, G minx2, G maxx2, G miny2,
                                 G maxy2, G minz2, G maxz2) {
                return !(maxx1 < minx2 || minx1 > maxx2 || maxy1 < miny2 || miny1 > maxy2 || maxz1 < minz2 ||
                         minz1 > maxz2);
            };

            // 3) Precompute local triangle AABBs (kept local, not communicated)
            std::vector<G> tminx(lnselements), tmaxx(lnselements), tminy(lnselements), tmaxy(lnselements),
                tminz(lnselements), tmaxz(lnselements);
#pragma omp parallel for
            for (ptrdiff_t i = 0; i < lnselements; ++i) {
                const I i0 = s0[i], i1 = s1[i], i2 = s2[i];
                const G x0 = sx[i0], x1 = sx[i1], x2 = sx[i2];
                const G y0 = sy[i0], y1 = sy[i1], y2 = sy[i2];
                const G z0 = sz[i0], z1 = sz[i1], z2 = sz[i2];
                tminx[i] = std::min({x0, x1, x2});
                tmaxx[i] = std::max({x0, x1, x2});
                tminy[i] = std::min({y0, y1, y2});
                tmaxy[i] = std::max({y0, y1, y2});
                tminz[i] = std::min({z0, z1, z2});
                tmaxz[i] = std::max({z0, z1, z2});
            }

            // 4) For each rank, compact intersecting triangles and their unique vertices
            std::vector<std::vector<G>> send_vx(comm_size), send_vy(comm_size), send_vz(comm_size);
            std::vector<std::vector<I>> send_s0(comm_size), send_s1(comm_size), send_s2(comm_size);
            std::vector<int> vmark(lnspoints, -1);
            std::vector<I> vmap_idx(lnspoints, 0);
            int mark_token = 0;

            for (int r = 0; r < comm_size; ++r) {
                const G bminx = all_pminx[r], bmaxx = all_pmaxx[r];
                const G bminy = all_pminy[r], bmaxy = all_pmaxy[r];
                const G bminz = all_pminz[r], bmaxz = all_pmaxz[r];

                // First pass: count intersecting triangles (no reallocs)
                ptrdiff_t tri_count = 0;
                for (ptrdiff_t tri = 0; tri < lnselements; ++tri) {
                    if (intersects(tminx[tri], tmaxx[tri], tminy[tri], tmaxy[tri], tminz[tri], tmaxz[tri], bminx, bmaxx,
                                   bminy, bmaxy, bminz, bmaxz)) {
                        ++tri_count;
                    }
                }

                printf("Rank %d found %ld/%ld intersecting triangles\n", comm_rank, tri_count, lnselements);
                if (tri_count == 0) continue;

                std::vector<I> tri_indices(static_cast<size_t>(tri_count));
                ptrdiff_t tri_cursor = 0;
                for (ptrdiff_t tri = 0; tri < lnselements; ++tri) {
                    if (intersects(tminx[tri], tmaxx[tri], tminy[tri], tmaxy[tri], tminz[tri], tmaxz[tri], bminx, bmaxx,
                                   bminy, bmaxy, bminz, bmaxz)) {
                        tri_indices[static_cast<size_t>(tri_cursor++)] = static_cast<I>(tri);
                    }
                }

                // Count unique vertices using reusable marker array
                ++mark_token;
                ptrdiff_t unique_count = 0;
                auto count_vertex = [&](I vid) {
                    if (vmark[vid] != mark_token) {
                        vmark[vid] = mark_token;
                        ++unique_count;
                    }
                };

                for (I tri_idx : tri_indices) {
                    count_vertex(s0[tri_idx]);
                    count_vertex(s1[tri_idx]);
                    count_vertex(s2[tri_idx]);
                }

                // Build vertex list and remap table
                std::vector<I> vlist(static_cast<size_t>(unique_count));

                // Second pass: assign ids in order of discovery
                ++mark_token;
                ptrdiff_t cursor = 0;
                auto add_vertex = [&](I vid) -> I {
                    if (vmark[vid] != mark_token) {
                        vmark[vid] = mark_token;
                        const I id = static_cast<I>(cursor);
                        vmap_idx[vid] = id;
                        vlist[static_cast<size_t>(cursor)] = vid;
                        ++cursor;
                        return id;
                    }
                    return vmap_idx[vid];
                };

                std::vector<I> tri_s0, tri_s1, tri_s2;
                tri_s0.reserve(tri_indices.size());
                tri_s1.reserve(tri_indices.size());
                tri_s2.reserve(tri_indices.size());

                for (I tri_idx : tri_indices) {
                    const I m0 = add_vertex(s0[tri_idx]);
                    const I m1 = add_vertex(s1[tri_idx]);
                    const I m2 = add_vertex(s2[tri_idx]);
                    tri_s0.push_back(m0);
                    tri_s1.push_back(m1);
                    tri_s2.push_back(m2);
                }

                auto &vx = send_vx[r];
                auto &vy = send_vy[r];
                auto &vz = send_vz[r];
                vx.resize(vlist.size());
                vy.resize(vlist.size());
                vz.resize(vlist.size());
                for (size_t idx = 0; idx < vlist.size(); ++idx) {
                    const I vid = vlist[idx];
                    vx[idx] = sx[vid];
                    vy[idx] = sy[vid];
                    vz[idx] = sz[vid];
                }

                auto &ts0 = send_s0[r];
                auto &ts1 = send_s1[r];
                auto &ts2 = send_s2[r];
                ts0 = std::move(tri_s0);
                ts1 = std::move(tri_s1);
                ts2 = std::move(tri_s2);
            }

            std::vector<int> send_vcounts(comm_size, 0), send_tcounts(comm_size, 0);
            for (int r = 0; r < comm_size; ++r) {
                send_vcounts[r] = static_cast<int>(send_vx[r].size());
                send_tcounts[r] = static_cast<int>(send_s0[r].size());
            }

            // 5) Exchange counts
            std::vector<int> recv_vcounts(comm_size, 0), recv_tcounts(comm_size, 0);
            MPI_Alltoall(send_vcounts.data(), 1, MPI_INT, recv_vcounts.data(), 1, MPI_INT, comm);
            MPI_Alltoall(send_tcounts.data(), 1, MPI_INT, recv_tcounts.data(), 1, MPI_INT, comm);

            // Build displacements
            std::vector<int> send_vdisp(comm_size, 0), recv_vdisp(comm_size, 0), send_tdisp(comm_size, 0),
                recv_tdisp(comm_size, 0);
            for (int r = 1; r < comm_size; ++r) {
                send_vdisp[r] = send_vdisp[r - 1] + send_vcounts[r - 1];
                recv_vdisp[r] = recv_vdisp[r - 1] + recv_vcounts[r - 1];
                send_tdisp[r] = send_tdisp[r - 1] + send_tcounts[r - 1];
                recv_tdisp[r] = recv_tdisp[r - 1] + recv_tcounts[r - 1];
            }

            const int total_recv_v = recv_vdisp.back() + recv_vcounts.back();
            const int total_recv_t = recv_tdisp.back() + recv_tcounts.back();

            // 6) Flatten send buffers for vertices (SoA) and triangles, then exchange
            const int total_send_v = send_vdisp.back() + send_vcounts.back();
            const int total_send_t = send_tdisp.back() + send_tcounts.back();

            std::vector<G> send_vbuf_x(total_send_v), send_vbuf_y(total_send_v), send_vbuf_z(total_send_v);
            std::vector<I> send_tbuf0(total_send_t), send_tbuf1(total_send_t), send_tbuf2(total_send_t);

            for (int r = 0; r < comm_size; ++r) {
                if (send_vcounts[r] > 0) {
                    std::copy(send_vx[r].begin(), send_vx[r].end(), send_vbuf_x.begin() + send_vdisp[r]);
                    std::copy(send_vy[r].begin(), send_vy[r].end(), send_vbuf_y.begin() + send_vdisp[r]);
                    std::copy(send_vz[r].begin(), send_vz[r].end(), send_vbuf_z.begin() + send_vdisp[r]);
                }
                if (send_tcounts[r] > 0) {
                    std::copy(send_s0[r].begin(), send_s0[r].end(), send_tbuf0.begin() + send_tdisp[r]);
                    std::copy(send_s1[r].begin(), send_s1[r].end(), send_tbuf1.begin() + send_tdisp[r]);
                    std::copy(send_s2[r].begin(), send_s2[r].end(), send_tbuf2.begin() + send_tdisp[r]);
                }
            }

            std::vector<G> recv_sx(total_recv_v), recv_sy(total_recv_v), recv_sz(total_recv_v);
            MPI_Alltoallv(send_vbuf_x.data(), send_vcounts.data(), send_vdisp.data(), mpi_type_G,
                          recv_sx.data(), recv_vcounts.data(), recv_vdisp.data(), mpi_type_G, comm);
            MPI_Alltoallv(send_vbuf_y.data(), send_vcounts.data(), send_vdisp.data(), mpi_type_G,
                          recv_sy.data(), recv_vcounts.data(), recv_vdisp.data(), mpi_type_G, comm);
            MPI_Alltoallv(send_vbuf_z.data(), send_vcounts.data(), send_vdisp.data(), mpi_type_G,
                          recv_sz.data(), recv_vcounts.data(), recv_vdisp.data(), mpi_type_G, comm);

            std::vector<I> recv_s0(total_recv_t), recv_s1(total_recv_t), recv_s2(total_recv_t);
            MPI_Alltoallv(send_tbuf0.data(), send_tcounts.data(), send_tdisp.data(), mpi_type_I,
                          recv_s0.data(), recv_tcounts.data(), recv_tdisp.data(), mpi_type_I, comm);
            MPI_Alltoallv(send_tbuf1.data(), send_tcounts.data(), send_tdisp.data(), mpi_type_I,
                          recv_s1.data(), recv_tcounts.data(), recv_tdisp.data(), mpi_type_I, comm);
            MPI_Alltoallv(send_tbuf2.data(), send_tcounts.data(), send_tdisp.data(), mpi_type_I,
                          recv_s2.data(), recv_tcounts.data(), recv_tdisp.data(), mpi_type_I, comm);

            // 8) Adjust indices in-place by per-source vertex offsets and run EDF on recv buffers
            for (int r = 0; r < comm_size; ++r) {
                const ptrdiff_t v_offset = static_cast<ptrdiff_t>(recv_vdisp[r]);
                const ptrdiff_t tri_offset = static_cast<ptrdiff_t>(recv_tdisp[r]);
                const ptrdiff_t tri_end = tri_offset + recv_tcounts[r];
                ptrdiff_t local_idx = recv_tdisp[r];
                for (ptrdiff_t t = tri_offset; t < tri_end; ++t, ++local_idx) {
                    recv_s0[local_idx] = static_cast<I>(recv_s0[local_idx] + v_offset);
                    recv_s1[local_idx] = static_cast<I>(recv_s1[local_idx] + v_offset);
                    recv_s2[local_idx] = static_cast<I>(recv_s2[local_idx] + v_offset);
                }
            }

            // if(!comm_rank) {
                printf("Rank %d received %lld triangles and %lld vertices\n", comm_rank, total_recv_t, total_recv_v);
            // }

            ssdf::edf_select<G, T, I>(lnpoints,
                               x,
                               y,
                               z,
                               static_cast<ptrdiff_t>(total_recv_t),
                               recv_s0.data(),
                               recv_s1.data(),
                               recv_s2.data(),
                               static_cast<ptrdiff_t>(total_recv_v),
                               recv_sx.data(),
                               recv_sy.data(),
                               recv_sz.data(),
                               out);

            return 0;
        
    }
}  // namespace ssdf