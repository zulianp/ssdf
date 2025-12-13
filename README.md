# Simple Signed Distance Function

## Building

### Prerequisites

- CMake 3.15 or higher
- C++17 compatible compiler

### Optional Dependencies

- **OpenMP**: For shared memory parallelization (Avaialble)
- **MPI**: For distributed memory parallelization (Available)
- **CUDA**: For GPU acceleration (TODO)
- **HIP**: For AMD GPU acceleration (TODO)

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
   
   Or with optional dependencies:
   ```bash
   cmake -DSSDF_ENABLE_MPI=ON -DSSDF_ENABLE_OPENMP=ON ..
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
- `SSDF_ENABLE_HIP`: Enable HIP support (default: OFF)

### Using ssdf in Your CMake Project

After installation, you can use ssdf in your CMake project:

```cmake
find_package(ssdf REQUIRED)
target_link_libraries(your_target PRIVATE ssdf::ssdf)
```


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

- API: `ssdf::edf<G,T,I>(npoints, x, y, z, nselements, s0, s1, s2, nspoints, sx, sy, sz, out)`.
- Inputs: point SoA `(x,y,z)`, surface vertex SoA `(sx,sy,sz)`, triangle indices `(s0,s1,s2)`.
- Output: `out[i]` is the (unsigned) distance from point i to the surface.
- Binary driver: `sdf_exe <surf_folder> <points_folder> <output_file>`
  - `surf_folder`: `x.raw,y.raw,z.raw` (float), `i0.raw,i1.raw,i2.raw` (int)
  - `points_folder`: `x.raw,y.raw,z.raw` (float)
  - Output: binary float array of distances

## Using the MPI SDF (pedf)

- API: `ssdf::pedf<G,T,I>(comm, lnpoints, x, y, z, lnselements, s0, s1, s2, lnspoints, sx, sy, sz, out, use_allgather=true)`.
- Data distribution: points and surface are partitioned across ranks.
- Binary driver: `psdf_exe <surf_folder> <points_folder> <output_file>`
  - Reads SoA `x.raw,y.raw,z.raw` and `i0.raw,i1.raw,i2.raw` via `MPI_File_read_at_all`
  - Writes distributed output with `MPI_File_write_at_all`
  - Optional `SSDF_INPUT_SDF` to initialize output;


## Usage within FSI simulations

In order to save computational time only compute for what matters

1) Compute the edf for the fluid domain keep the `edf` untouched for the rest of the simulation
2) Create a new array copy of `edf` update it with the current deformed structural surface