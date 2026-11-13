# Simple Signed Distance Function

## Building

### Prerequisites

- CMake 3.15 or higher
- C++17 compatible compiler

### Optional Dependencies

- **OpenMP**: For shared memory parallelization (Avaialble)
- **MPI**: For distributed memory parallelization (Available)
- **CUDA**: For GPU acceleration (Rough integration with cuBQL)
- **HIP**: For AMD GPUs: the same device BVH, compiled with HIP (see [AMD GPUs (HIP)](#amd-gpus-hip))
- **cuBQL**: For BVH datastructures, fastest implementation a.t.m. (https://github.com/NVIDIA/cuBQL.git)

### Building the Project

1. **Clone the repository:**
   ```bash
   git clone https://github.com/zulianp/ssdf.git
   cd ssdf
   ```

2. **Create a build directory:**
   ```bash
   mkdir build
   cd build
   ```

3. **Configure with CMake:**
   ```bash
   cmake ..
   ```
   
   With optional dependencies:
   ```bash
   cmake -DSSDF_ENABLE_MPI=ON -DSSDF_ENABLE_OPENMP=ON ..
   ```

   With CUBQL
   ```bash
   cmake -DSSDF_ENABLE_MPI=ON -DSSDF_ENABLE_OPENMP=ON -DSSDF_ENABLE_CUBIQL=ON -DcuBQL_DIR=<CUBQL_PROJECT_FOLDER> ..
   ```


4. **Build:**
   ```bash
   cmake --build .
   ```

5. **Install (optional):**
   ```bash
   cmake --install . --prefix /path/to/install
   ```

### CMake Options

- `BUILD_SHARED_LIBS`: Build shared library (default: ON)
- `SSDF_ENABLE_MPI`: Enable MPI support (default: OFF)
- `SSDF_ENABLE_OPENMP`: Enable OpenMP support (default: OFF)
- `SSDF_ENABLE_CUDA`: Enable CUDA support (default: OFF)
- `SSDF_ENABLE_HIP`: With `SSDF_ENABLE_CUDA=ON`, compile the device BVH with HIP instead of CUDA, for AMD GPUs (default: OFF)
- `SSDF_HIP_ARCHITECTURES`: AMD GPU architectures, used when `CMAKE_HIP_ARCHITECTURES` is not set (default: gfx942)
- `SSDF_ENABLE_CUBIQL`: Use cuBQL for the BVH (default: ON)

### AMD GPUs (HIP)

The device BVH (`serial/bvh/bvh.cu`, built on cuBQL) is written against the CUDA runtime API and
compiles with HIP. `SSDF_ENABLE_CUDA=ON` turns the device path on; adding `SSDF_ENABLE_HIP=ON`
compiles it with HIP instead of CUDA, together with cuBQL's GPU builder, which sorts with hipCUB.

Requirements: ROCm with HIP, hipCUB and rocPRIM, CMake 3.21 or newer (the first version with the
HIP language), and the cuBQL submodule (`git submodule update --init external/cubql`).

```bash
cmake -S . -B build-hip -DCMAKE_BUILD_TYPE=Release \
    -DSSDF_ENABLE_CUDA=ON -DSSDF_ENABLE_HIP=ON -DSSDF_ENABLE_CUBIQL=ON \
    -DSSDF_HIP_ARCHITECTURES=gfx942 \
    -DCMAKE_HIP_COMPILER=/path/to/rocm/clang++
cmake --build build-hip -j 8
```

- **Architectures:** set `SSDF_HIP_ARCHITECTURES` to your GPU's architecture, for example `gfx942`
  for MI300A/MI300X or `gfx90a` for MI250X; `rocminfo | grep gfx` prints it. A parent project's
  `CMAKE_HIP_ARCHITECTURES` takes precedence.
- **Compilers:** only the `.cu` files are compiled by the HIP compiler, ROCm's `clang++`. In a
  standard ROCm install it is `$ROCM_PATH/llvm/bin/clang++`; module systems often put it
  elsewhere, so check with `which clang++` after loading ROCm. The other sources are compiled by
  the C++ compiler (for example GCC), and the library links `hip::host` publicly.
- **Device memory:** `serial/bvh/ssdf_gpu_runtime.cuh` maps the CUDA names that `bvh.cu` and cuBQL
  call onto HIP. It declares the CUDA 12 runtime API, so cuBQL allocates with `hipMallocAsync`
  instead of managed memory, which an AMD GPU running without XNACK cannot page in.
- **cuBQL:** `CUBQL_DISABLE_CUDA` is set, so cuBQL does not enable CUDA on its own when it finds
  nvcc. Its headers are compiled by ssdf's HIP sources.
- **Floating point:** subnormals are kept (`-fno-gpu-flush-denormals-to-zero`), as nvcc keeps them
  by default.
- **As a subproject:** set `SSDF_ENABLE_CUDA` and `SSDF_ENABLE_HIP` before `add_subdirectory`. If
  the parent has already enabled the HIP language, ssdf uses its settings. SFEM does this with
  `SFEM_ENABLE_HIP`.

### Using ssdf in Your CMake Project

After installation, you can use ssdf in your CMake project:

```cmake
find_package(ssdf REQUIRED)
target_link_libraries(your_target PRIVATE ssdf::ssdf)
```

The build also generates a flat, installable configuration header:

- `#include <ssdf_config.hpp>`: compile-time feature macros (MPI/OpenMP/CUDA/...)


# Cite SSDF

Cite SSDF if you use it for your work:

```bibtex
@misc{SSDFgit,
	author = {Zulian, Patrick,
	title = {{SSDF}: Simple {Signed Distance Fields}},
	url = {https://github.com/zulianp/ssdf},
	howpublished = {https://github.com/zulianp/ssdf},
	year = {2025}
}
```

## Using the CPU SDF (edf)

- API: `ssdf::edf<G,T,I>(...)` (reference), `ssdf::edf_celllist<G,T,I>(...)` (cell-list, best performance version)
- Use `ssdf::edf_select<G,T,I>(...)` for using the best perfoming version selected by the implementor
- Inputs: point SoA `(x,y,z)`, surface vertex SoA `(sx,sy,sz)`, triangle indices `(s0,s1,s2)`.
- Output: `out[i]` is the (unsigned) distance from point i to the surface.
- Cell-list variant is untesed for now use edf for correctness
- Optional validation: set `SSDF_CELL_VALIDATE=1` alongside to run both paths and report max difference (copies the cell-list result to `out`).
- Binary driver: `sdf <surf_folder> <points_folder> <output_file>`
  - `surf_folder`: `x.raw,y.raw,z.raw` (float), `i0.raw,i1.raw,i2.raw` (int)
  - `points_folder`: `x.raw,y.raw,z.raw` (float)
  - Output: binary float array of distances

## Using the MPI SDF (pedf)

- API: `ssdf::pedf<G,T,I>(comm, lnpoints, x, y, z, lnselements, s0, s1, s2, lnspoints, sx, sy, sz, out)`.
- Data distribution: points and surface are partitioned across ranks.
- Binary driver: `psdf <surf_folder> <points_folder> <output_file>` (requires `-DSSDF_ENABLE_MPI=ON`)
  - Reads SoA `x.raw,y.raw,z.raw` and `i0.raw,i1.raw,i2.raw` via `MPI_File_read_at_all`
  - Writes distributed output with `MPI_File_write_at_all`
  - Optional `SSDF_INPUT_SDF` to initialize output;

## Using the CPU squared-distance SDF (sqedf)

- API: `ssdf::sqedf<G,T,I>(...)` (works on squared distances internally).
- Binary driver: `sqedf <surf_folder1> [surf_folder2] ... <points_folder> <output_file>`


## Usage within FSI simulations

In order to save computational time only compute for what matters

1) Compute the edf for the fluid domain keep the `edf` untouched for the rest of the simulation
2) Create a new array copy of `edf` update it with the current deformed structural surface