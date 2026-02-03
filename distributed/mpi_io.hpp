#ifndef SSDF_DISTRIBUTED_MPI_IO_HPP
#define SSDF_DISTRIBUTED_MPI_IO_HPP

#include "psdf.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

namespace ssdf {
namespace mpi_io {

// Collective distributed read: each rank gets a contiguous chunk.
// - global_count_out: total number of elements in the file
// - offset_out: element offset of this rank's chunk
template <typename T>
inline ptrdiff_t read_distributed(const std::string &path,
                                  std::vector<T> &local,
                                  MPI_Comm comm,
                                  ptrdiff_t &global_count_out,
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
    const ptrdiff_t elem_size = static_cast<ptrdiff_t>(sizeof(T));
    if (fsize % elem_size != 0) {
        if (rank == 0) std::fprintf(stderr, "Error: size mismatch for %s\n", path.c_str());
        MPI_Abort(comm, 1);
    }
    global_count_out = static_cast<ptrdiff_t>(fsize / elem_size);

    const ptrdiff_t base = global_count_out / size;
    const ptrdiff_t rem = global_count_out % size;
    const ptrdiff_t local_count = base + (rank < rem ? 1 : 0);
    const ptrdiff_t offset = base * rank + std::min<ptrdiff_t>(rank, rem);
    offset_out = offset;

    local.resize(static_cast<size_t>(local_count));
    const MPI_Offset disp = static_cast<MPI_Offset>(offset * elem_size);

    const MPI_Datatype dtype = ssdf::mpi_type<T>();
    if (local_count > static_cast<ptrdiff_t>(std::numeric_limits<int>::max())) {
        if (rank == 0) std::fprintf(stderr, "Error: local_count overflows MPI int count for %s\n", path.c_str());
        MPI_Abort(comm, 1);
    }

    MPI_File_set_view(fh, disp, dtype, dtype, "native", MPI_INFO_NULL);
    MPI_Status st;
    MPI_File_read_at_all(fh, 0, local.data(), static_cast<int>(local_count), dtype, &st);
    MPI_File_close(&fh);
    return local_count;
}

// Write distributed array (contiguous partition) to file.
template <typename T>
inline void write_distributed(const std::string &path, const std::vector<T> &local, ptrdiff_t offset, MPI_Comm comm) {
    int rank = 0;
    MPI_Comm_rank(comm, &rank);

    MPI_File fh;
    if (MPI_File_open(comm, path.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh) != MPI_SUCCESS) {
        if (rank == 0) std::fprintf(stderr, "Error: cannot open output %s\n", path.c_str());
        MPI_Abort(comm, 1);
    }

    const MPI_Datatype dtype = ssdf::mpi_type<T>();
    if (local.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
        if (rank == 0) std::fprintf(stderr, "Error: local.size() overflows MPI int count for %s\n", path.c_str());
        MPI_Abort(comm, 1);
    }

    const MPI_Offset disp = static_cast<MPI_Offset>(offset * static_cast<ptrdiff_t>(sizeof(T)));
    MPI_File_set_view(fh, disp, dtype, dtype, "native", MPI_INFO_NULL);
    MPI_Status st;
    MPI_File_write_at_all(fh, 0, local.data(), static_cast<int>(local.size()), dtype, &st);
    MPI_File_close(&fh);
}

}  // namespace mpi_io
}  // namespace ssdf

#endif  // SSDF_DISTRIBUTED_MPI_IO_HPP
