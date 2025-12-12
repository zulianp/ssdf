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
    const MPI_Datatype dtype = mpi_type<T>();
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
void mpi_write_distributed(const std::string &path,
                           const std::vector<T> &local,
                           ptrdiff_t offset,
                           MPI_Comm comm) {
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    MPI_File fh;
    if (MPI_File_open(comm, path.c_str(),
                      MPI_MODE_CREATE | MPI_MODE_WRONLY,
                      MPI_INFO_NULL, &fh) != MPI_SUCCESS) {
        if (rank == 0) std::fprintf(stderr, "Error: cannot open output %s\n", path.c_str());
        MPI_Abort(comm, 1);
    }
    const MPI_Datatype dtype = mpi_type<T>();
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
            std::fprintf(stderr,
                         "Usage: %s <surf_folder> <points_folder> <output_file>\n",
                         argv[0]);
            std::fprintf(stderr, "Each folder must contain x.raw, y.raw, z.raw (float) and i0.raw, i1.raw, i2.raw (int) for surface.\n");
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
        if (!folder.empty() && folder.back() == '/')
            return folder + fname;
        return folder + "/" + fname;
    };

    // Read surface points (distributed)
    ptrdiff_t surf_nspoints = 0, surf_sp_offset = 0;
    std::vector<G> sx, sy, sz;
    mpi_read_distributed<G>(make_path(surf_folder, "x.raw"), sx, comm, surf_nspoints, surf_sp_offset);
    ptrdiff_t tmp_global = 0, tmp_off = 0;
    mpi_read_distributed<G>(make_path(surf_folder, "y.raw"), sy, comm, tmp_global, tmp_off);
    mpi_read_distributed<G>(make_path(surf_folder, "z.raw"), sz, comm, tmp_global, tmp_off);

    // Read surface indices (distributed), adjust to local indexing
    ptrdiff_t surf_nselements = 0, surf_se_offset = 0;
    std::vector<I> s0, s1, s2;
    mpi_read_distributed<I>(make_path(surf_folder, "i0.raw"), s0, comm, surf_nselements, surf_se_offset);
    mpi_read_distributed<I>(make_path(surf_folder, "i1.raw"), s1, comm, tmp_global, tmp_off);
    mpi_read_distributed<I>(make_path(surf_folder, "i2.raw"), s2, comm, tmp_global, tmp_off);

    // Convert global indices to local by subtracting surf_sp_offset
    for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(s0.size()); ++i) {
        s0[i] -= static_cast<I>(surf_sp_offset);
        s1[i] -= static_cast<I>(surf_sp_offset);
        s2[i] -= static_cast<I>(surf_sp_offset);
    }

    // Read points (distributed)
    ptrdiff_t npoints_global = 0, points_offset = 0;
    std::vector<G> x, y, z;
    mpi_read_distributed<G>(make_path(points_folder, "x.raw"), x, comm, npoints_global, points_offset);
    mpi_read_distributed<G>(make_path(points_folder, "y.raw"), y, comm, tmp_global, tmp_off);
    mpi_read_distributed<G>(make_path(points_folder, "z.raw"), z, comm, tmp_global, tmp_off);
    const ptrdiff_t lnpoints = static_cast<ptrdiff_t>(x.size());

    // Output buffer
    

    // Choose mode via env PEDF_USE_ALLGATHER (default true)
    bool use_allgather = true;
    if (const char *env = std::getenv("PEDF_USE_ALLGATHER")) {
        use_allgather = (std::strcmp(env, "0") != 0 && std::strcmp(env, "false") != 0 && std::strcmp(env, "FALSE") != 0);
    }

    const char *INPUT_SDF = std::getenv("INPUT_SDF");

    std::vector<T> out;
    if (INPUT_SDF) { 
        mpi_read_distributed<T>(INPUT_SDF, out, comm, npoints_global, points_offset);
    }
    else {
        out.resize(lnpoints, std::numeric_limits<T>::max());
    }
    // Compute EDF for local points vs distributed surface
    ssdf::pedf<G, T, I>(comm,
                        lnpoints,
                        x.data(), y.data(), z.data(),
                        static_cast<ptrdiff_t>(s0.size()), s0.data(), s1.data(), s2.data(),
                        static_cast<ptrdiff_t>(sx.size()), sx.data(), sy.data(), sz.data(),
                        out.data(),
                        use_allgather);

    // Write output
    mpi_write_distributed<T>(output_file, out, points_offset, comm);

    MPI_Finalize();
    return 0;
}