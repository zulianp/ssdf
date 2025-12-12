# Simple Signed Distance Function

## Building

### Prerequisites

- CMake 3.15 or higher
- C++17 compatible compiler

### Optional Dependencies

- **OpenMP**: For shared memory parallelization
- **MPI**: For distributed memory parallelization (TODO
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
