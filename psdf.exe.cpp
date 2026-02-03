// Parallel SDF driver using MPI I/O and pedf()
// - Points are distributed across ranks via MPI_File_read_at_all
// - Surface mesh is distributed across ranks (contiguous partition)
// - pedf() computes SDF for local points
// - Output is written with MPI_File_write_at_all (contiguous partition)

#include <mpi.h>
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "psdf.hpp"

namespace {

    // Collective distributed read: each rank gets a contiguous chunk
    template <typename T>
    ptrdiff_t mpi_read_distributed(const std::string &path,
                                   std::vector<T> &local,
                                   MPI_Comm comm,
                                   ptrdiff_t &global_count,
                                   ptrdiff_t &offset_out) {
        int rank = 0, size = 0;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &size);

        MPI_File fh;
        if (MPI_File_open(comm, path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh) != MPI_SUCCESS) {
            if (rank == 0) std::fprintf(stderr, "Error: cannot open %s\n", path.c_str());
            MPI_Abort(comm, 1);
        }

        MPI_Offset fsize = 0;
        MPI_File_get_size(fh, &fsize);
        const MPI_Datatype dtype = ssdf::mpi_type<T>();
        const ptrdiff_t elem_size = static_cast<ptrdiff_t>(sizeof(T));
        if (fsize % elem_size != 0) {
            if (rank == 0) std::fprintf(stderr, "Error: size mismatch for %s\n", path.c_str());
            MPI_Abort(comm, 1);
        }
        global_count = static_cast<ptrdiff_t>(fsize / elem_size);

        const ptrdiff_t base = global_count / size;
        const ptrdiff_t rem = global_count % size;
        const ptrdiff_t local_count = base + (rank < rem ? 1 : 0);
        const ptrdiff_t offset = base * rank + std::min<ptrdiff_t>(rank, rem);
        offset_out = offset;

        local.resize(local_count);
        const MPI_Offset disp = static_cast<MPI_Offset>(offset * elem_size);
        MPI_File_set_view(fh, disp, dtype, dtype, "native", MPI_INFO_NULL);
        MPI_Status st;
        MPI_File_read_at_all(fh, 0, local.data(), static_cast<int>(local_count), dtype, &st);
        MPI_File_close(&fh);
        return local_count;
    }

    // Write distributed array (contiguous partition) to file
    template <typename T>
    void mpi_write_distributed(const std::string &path, const std::vector<T> &local, ptrdiff_t offset, MPI_Comm comm) {
        int rank = 0;
        MPI_Comm_rank(comm, &rank);
        MPI_File fh;
        if (MPI_File_open(comm, path.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh) != MPI_SUCCESS) {
            if (rank == 0) std::fprintf(stderr, "Error: cannot open output %s\n", path.c_str());
            MPI_Abort(comm, 1);
        }
        const MPI_Datatype dtype = ssdf::mpi_type<T>();
        const MPI_Offset disp = static_cast<MPI_Offset>(offset * static_cast<ptrdiff_t>(sizeof(T)));
        MPI_File_set_view(fh, disp, dtype, dtype, "native", MPI_INFO_NULL);
        MPI_Status st;
        MPI_File_write_at_all(fh, 0, local.data(), static_cast<int>(local.size()), dtype, &st);
        MPI_File_close(&fh);
    }

}  // namespace

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank = 0;
    MPI_Comm_rank(comm, &rank);

    if (argc != 4) {
        if (rank == 0) {
            std::fprintf(stderr, "Usage: %s <surf_folder> <points_folder> <output_file>\n", argv[0]);
            std::fprintf(
                stderr,
                "Each folder must contain x.raw, y.raw, z.raw (float) and i0.raw, i1.raw, i2.raw (int) for surface.\n");
        }
        MPI_Finalize();
        return 1;
    }

    const std::string surf_folder = argv[1];
    const std::string points_folder = argv[2];
    const std::string output_file = argv[3];

    using G = float;
    using T = float;
    using I = int;

    auto make_path = [](const std::string &folder, const char *fname) {
        if (!folder.empty() && folder.back() == '/') return folder + fname;
        return folder + "/" + fname;
    };

    if (!rank) {
        printf("Reading mesh files!\n");
        fflush(stdout);
    }

    // Read surface points (distributed)
    ptrdiff_t surf_nspoints = 0, surf_sp_offset = 0;
    std::vector<G> sx, sy, sz;
    mpi_read_distributed<G>(make_path(surf_folder, "x.raw"), sx, comm, surf_nspoints, surf_sp_offset);
    ptrdiff_t tmp_global = 0, tmp_off = 0;
    mpi_read_distributed<G>(make_path(surf_folder, "y.raw"), sy, comm, tmp_global, tmp_off);
    mpi_read_distributed<G>(make_path(surf_folder, "z.raw"), sz, comm, tmp_global, tmp_off);

    // Read surface indices (distributed)
    ptrdiff_t surf_nselements = 0, surf_se_offset = 0;
    std::vector<I> s0, s1, s2;
    mpi_read_distributed<I>(make_path(surf_folder, "i0.raw"), s0, comm, surf_nselements, surf_se_offset);
    mpi_read_distributed<I>(make_path(surf_folder, "i1.raw"), s1, comm, tmp_global, tmp_off);
    mpi_read_distributed<I>(make_path(surf_folder, "i2.raw"), s2, comm, tmp_global, tmp_off);

    // Compact indices and points: gather all surface points, then keep only those referenced by local triangles
    {
        int size = 0;
        MPI_Comm_size(comm, &size);

        // Gather all surface points to all ranks
        std::vector<long long> all_spoints(size);
        long long lnspoints_ll = static_cast<long long>(sx.size());
        MPI_Allgather(&lnspoints_ll, 1, MPI_LONG_LONG, all_spoints.data(), 1, MPI_LONG_LONG, comm);

        std::vector<long long> off_spoints(size + 1, 0);
        for (int r = 0; r < size; ++r) {
            off_spoints[r + 1] = off_spoints[r] + all_spoints[r];
        }
        const ptrdiff_t nspoints_total = static_cast<ptrdiff_t>(off_spoints.back());

        std::vector<G> gx(nspoints_total), gy(nspoints_total), gz(nspoints_total);
        std::vector<int> sc_sp(size), dsp_sp(size);
        for (int r = 0; r < size; ++r) {
            sc_sp[r] = static_cast<int>(all_spoints[r]);
            dsp_sp[r] = static_cast<int>(off_spoints[r]);
        }

        const MPI_Datatype dtype_G = ssdf::mpi_type<G>();
        MPI_Allgatherv(sx.data(), sc_sp[rank], dtype_G, gx.data(), sc_sp.data(), dsp_sp.data(), dtype_G, comm);
        MPI_Allgatherv(sy.data(), sc_sp[rank], dtype_G, gy.data(), sc_sp.data(), dsp_sp.data(), dtype_G, comm);
        MPI_Allgatherv(sz.data(), sc_sp[rank], dtype_G, gz.data(), sc_sp.data(), dsp_sp.data(), dtype_G, comm);

        // Find unique point indices referenced by local triangles
        std::vector<bool> point_used(nspoints_total, false);
        for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(s0.size()); ++i) {
            if (s0[i] >= 0 && s0[i] < static_cast<I>(nspoints_total)) {
                point_used[s0[i]] = true;
            }
            if (s1[i] >= 0 && s1[i] < static_cast<I>(nspoints_total)) {
                point_used[s1[i]] = true;
            }
            if (s2[i] >= 0 && s2[i] < static_cast<I>(nspoints_total)) {
                point_used[s2[i]] = true;
            }
        }

        // Create mapping from global index to compacted local index
        std::vector<I> global_to_local(nspoints_total, -1);
        ptrdiff_t compacted_count = 0;
        for (ptrdiff_t i = 0; i < nspoints_total; ++i) {
            if (point_used[i]) {
                global_to_local[i] = static_cast<I>(compacted_count);
                ++compacted_count;
            }
        }

        // Compact surface points
        sx.resize(compacted_count);
        sy.resize(compacted_count);
        sz.resize(compacted_count);
        for (ptrdiff_t i = 0; i < nspoints_total; ++i) {
            if (point_used[i]) {
                const ptrdiff_t local_idx = global_to_local[i];
                sx[local_idx] = gx[i];
                sy[local_idx] = gy[i];
                sz[local_idx] = gz[i];
            }
        }

        // Remap triangle indices to compacted local indices
        for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(s0.size()); ++i) {
            if (s0[i] >= 0 && s0[i] < static_cast<I>(nspoints_total)) {
                s0[i] = global_to_local[s0[i]];
            }
            if (s1[i] >= 0 && s1[i] < static_cast<I>(nspoints_total)) {
                s1[i] = global_to_local[s1[i]];
            }
            if (s2[i] >= 0 && s2[i] < static_cast<I>(nspoints_total)) {
                s2[i] = global_to_local[s2[i]];
            }
        }
    }

    // Read points (distributed)
    ptrdiff_t npoints_global = 0, points_offset = 0;
    std::vector<G> x, y, z;
    mpi_read_distributed<G>(make_path(points_folder, "x.raw"), x, comm, npoints_global, points_offset);
    mpi_read_distributed<G>(make_path(points_folder, "y.raw"), y, comm, tmp_global, tmp_off);
    mpi_read_distributed<G>(make_path(points_folder, "z.raw"), z, comm, tmp_global, tmp_off);
    const ptrdiff_t lnpoints = static_cast<ptrdiff_t>(x.size());

    // Output buffer

    bool isupdate = false;

    const char *SSDF_INPUT_SDF = std::getenv("SSDF_INPUT_SDF");

    std::vector<T> out;
    if (SSDF_INPUT_SDF) {
        mpi_read_distributed<T>(SSDF_INPUT_SDF, out, comm, npoints_global, points_offset);
        isupdate = true;
    } else {
        int SSDF_INIT_WITH_VOLUME_BOX = 0;
        SSDF_READ_ENV(SSDF_INIT_WITH_VOLUME_BOX, std::stoi);
        if (SSDF_INIT_WITH_VOLUME_BOX) {
            out.resize(lnpoints);
            G xmin, xmax, ymin, ymax, zmin, zmax;
            ssdf::compute_aabb(comm, lnpoints, x.data(), y.data(), z.data(), &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);
            ssdf::all_points_aabb_signed_distance<G, T>(
                lnpoints, x.data(), y.data(), z.data(), xmin, xmax, ymin, ymax, zmin, zmax, out.data());

#pragma omp parallel for
            for (ptrdiff_t i = 0; i < lnpoints; i++) {
                out[i] = std::abs(out[i]);
            }

            isupdate = true;
        } else {
            out.resize(lnpoints, std::numeric_limits<T>::max());
        }
    }

    if (!rank) {
        printf("Running pedf with %ld points and %ld surface elements\n", npoints_global, surf_nselements);
    }

    int SSDF_USE_EDF_UPDATE = 0;
    SSDF_READ_ENV(SSDF_USE_EDF_UPDATE, std::stoi);
    isupdate = SSDF_USE_EDF_UPDATE && isupdate;

    int SSDF_DOUBLE_PRECISION = 0;
    SSDF_READ_ENV(SSDF_DOUBLE_PRECISION, std::stoi);
    if (SSDF_DOUBLE_PRECISION == 1 && !std::is_same<T, double>::value) {
        // Convert input and out to double
        std::vector<double> x_double(x.size());
        std::vector<double> y_double(y.size());
        std::vector<double> z_double(z.size());
        std::vector<double> out_double(out.size());
        std::vector<double> sx_double(sx.size());
        std::vector<double> sy_double(sy.size());
        std::vector<double> sz_double(sz.size());

        for (ptrdiff_t i = 0; i < lnpoints; ++i) {
            x_double[i] = double(x[i]);
            y_double[i] = double(y[i]);
            z_double[i] = double(z[i]);
            out_double[i] = double(out[i]);
        }

        for (ptrdiff_t i = 0; i < sx.size(); ++i) {
            sx_double[i] = double(sx[i]);
            sy_double[i] = double(sy[i]);
            sz_double[i] = double(sz[i]);
        }

        if (isupdate) {
            ssdf::pedf_update<double, double, I>(comm,
                                                 lnpoints,
                                                 x_double.data(),
                                                 y_double.data(),
                                                 z_double.data(),
                                                 static_cast<ptrdiff_t>(s0.size()),
                                                 s0.data(),
                                                 s1.data(),
                                                 s2.data(),
                                                 static_cast<ptrdiff_t>(sx.size()),
                                                 sx_double.data(),
                                                 sy_double.data(),
                                                 sz_double.data(),
                                                 out_double.data());
        } else {
            ssdf::pedf<double, double, I>(comm,
                                          lnpoints,
                                          x_double.data(),
                                          y_double.data(),
                                          z_double.data(),
                                          static_cast<ptrdiff_t>(s0.size()),
                                          s0.data(),
                                          s1.data(),
                                          s2.data(),
                                          static_cast<ptrdiff_t>(sx.size()),
                                          sx_double.data(),
                                          sy_double.data(),
                                          sz_double.data(),
                                          out_double.data());
        }

        // Convert output to float
        for (ptrdiff_t i = 0; i < lnpoints; ++i) {
            out[i] = static_cast<T>(out_double[i]);
        }

    } else {
        if (isupdate) {
            ssdf::pedf_update<G, T, I>(comm,
                                       lnpoints,
                                       x.data(),
                                       y.data(),
                                       z.data(),
                                       static_cast<ptrdiff_t>(s0.size()),
                                       s0.data(),
                                       s1.data(),
                                       s2.data(),
                                       static_cast<ptrdiff_t>(sx.size()),
                                       sx.data(),
                                       sy.data(),
                                       sz.data(),
                                       out.data());
        } else {
            ssdf::pedf<G, T, I>(comm,
                                lnpoints,
                                x.data(),
                                y.data(),
                                z.data(),
                                static_cast<ptrdiff_t>(s0.size()),
                                s0.data(),
                                s1.data(),
                                s2.data(),
                                static_cast<ptrdiff_t>(sx.size()),
                                sx.data(),
                                sy.data(),
                                sz.data(),
                                out.data());
        }
    }

    // Write output
    mpi_write_distributed<T>(output_file, out, points_offset, comm);

    MPI_Finalize();
    return 0;
}