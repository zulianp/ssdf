#include "ssdf.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

// Helper function to get file size in bytes
static long get_file_size(FILE* f) {
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    fseek(f, 0, SEEK_SET);
    return size;
}

// Helper function to read raw file and return size
static ptrdiff_t read_raw_file(const char* filepath, void* data, size_t element_size) {
    FILE* f = fopen(filepath, "rb");
    if (!f) {
        fprintf(stderr, "Error: Cannot open file '%s'\n", filepath);
        return -1;
    }
    
    long file_size = get_file_size(f);
    if (file_size < 0) {
        fprintf(stderr, "Error: Cannot determine size of file '%s'\n", filepath);
        fclose(f);
        return -1;
    }
    
    if (file_size % element_size != 0) {
        fprintf(stderr, "Error: File '%s' size (%ld) is not a multiple of element size (%zu)\n", 
                filepath, file_size, element_size);
        fclose(f);
        return -1;
    }
    
    ptrdiff_t count = file_size / element_size;
    if (fread(data, element_size, count, f) != count) {
        fprintf(stderr, "Error: Failed to read file '%s'\n", filepath);
        fclose(f);
        return -1;
    }
    
    fclose(f);
    return count;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <surf_folder> <points_folder> <output_file>\n", argv[0]);
        fprintf(stderr, "\n");
        fprintf(stderr, "Surface folder should contain:\n");
        fprintf(stderr, "  - x.raw (float array)\n");
        fprintf(stderr, "  - y.raw (float array)\n");
        fprintf(stderr, "  - z.raw (float array)\n");
        fprintf(stderr, "  - i0.raw (int array)\n");
        fprintf(stderr, "  - i1.raw (int array)\n");
        fprintf(stderr, "  - i2.raw (int array)\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Points folder should contain:\n");
        fprintf(stderr, "  - x.raw (float array)\n");
        fprintf(stderr, "  - y.raw (float array)\n");
        fprintf(stderr, "  - z.raw (float array)\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Output file format:\n");
        fprintf(stderr, "  - float sdf[npoints] (binary array of SDF values)\n");
        return 1;
    }

    const char* surf_folder = argv[1];
    const char* points_folder = argv[2];
    const char* output_file = argv[3];

    using G = float;
    using T = float;
    using I = int;

    // Helper function to construct file paths
    auto make_path = [](const char* folder, const char* filename) -> std::string {
        std::string path = folder;
        if (path.back() != '/') {
            path += '/';
        }
        path += filename;
        return path;
    };

    // Read surface point coordinates
    std::string surf_x_path = make_path(surf_folder, "x.raw");
    std::string surf_y_path = make_path(surf_folder, "y.raw");
    std::string surf_z_path = make_path(surf_folder, "z.raw");
    
    // First, determine size by reading one file
    FILE* test_file = fopen(surf_x_path.c_str(), "rb");
    if (!test_file) {
        fprintf(stderr, "Error: Cannot open surface file '%s'\n", surf_x_path.c_str());
        return 1;
    }
    long file_size = get_file_size(test_file);
    fclose(test_file);
    
    if (file_size < 0 || file_size % sizeof(G) != 0) {
        fprintf(stderr, "Error: Invalid surface x.raw file size\n");
        return 1;
    }
    
    ptrdiff_t nspoints = file_size / sizeof(G);
    
    G* sx = new G[nspoints];
    G* sy = new G[nspoints];
    G* sz = new G[nspoints];
    
    if (read_raw_file(surf_x_path.c_str(), sx, sizeof(G)) != nspoints ||
        read_raw_file(surf_y_path.c_str(), sy, sizeof(G)) != nspoints ||
        read_raw_file(surf_z_path.c_str(), sz, sizeof(G)) != nspoints) {
        fprintf(stderr, "Error: Failed to read surface coordinates\n");
        delete[] sx;
        delete[] sy;
        delete[] sz;
        return 1;
    }

    // Read surface element indices
    std::string surf_i0_path = make_path(surf_folder, "i0.raw");
    std::string surf_i1_path = make_path(surf_folder, "i1.raw");
    std::string surf_i2_path = make_path(surf_folder, "i2.raw");
    
    test_file = fopen(surf_i0_path.c_str(), "rb");
    if (!test_file) {
        fprintf(stderr, "Error: Cannot open surface index file '%s'\n", surf_i0_path.c_str());
        delete[] sx;
        delete[] sy;
        delete[] sz;
        return 1;
    }
    file_size = get_file_size(test_file);
    fclose(test_file);
    
    if (file_size < 0 || file_size % sizeof(I) != 0) {
        fprintf(stderr, "Error: Invalid surface i0.raw file size\n");
        delete[] sx;
        delete[] sy;
        delete[] sz;
        return 1;
    }
    
    ptrdiff_t nselements = file_size / sizeof(I);
    
    I* s0 = new I[nselements];
    I* s1 = new I[nselements];
    I* s2 = new I[nselements];
    
    if (read_raw_file(surf_i0_path.c_str(), s0, sizeof(I)) != nselements ||
        read_raw_file(surf_i1_path.c_str(), s1, sizeof(I)) != nselements ||
        read_raw_file(surf_i2_path.c_str(), s2, sizeof(I)) != nselements) {
        fprintf(stderr, "Error: Failed to read surface indices\n");
        delete[] sx;
        delete[] sy;
        delete[] sz;
        delete[] s0;
        delete[] s1;
        delete[] s2;
        return 1;
    }

    // Read point coordinates
    std::string points_x_path = make_path(points_folder, "x.raw");
    std::string points_y_path = make_path(points_folder, "y.raw");
    std::string points_z_path = make_path(points_folder, "z.raw");
    
    test_file = fopen(points_x_path.c_str(), "rb");
    if (!test_file) {
        fprintf(stderr, "Error: Cannot open points file '%s'\n", points_x_path.c_str());
        delete[] sx;
        delete[] sy;
        delete[] sz;
        delete[] s0;
        delete[] s1;
        delete[] s2;
        return 1;
    }
    file_size = get_file_size(test_file);
    fclose(test_file);
    
    if (file_size < 0 || file_size % sizeof(G) != 0) {
        fprintf(stderr, "Error: Invalid points x.raw file size\n");
        delete[] sx;
        delete[] sy;
        delete[] sz;
        delete[] s0;
        delete[] s1;
        delete[] s2;
        return 1;
    }
    
    ptrdiff_t npoints = file_size / sizeof(G);
    
    G* x = new G[npoints];
    G* y = new G[npoints];
    G* z = new G[npoints];
    
    if (read_raw_file(points_x_path.c_str(), x, sizeof(G)) != npoints ||
        read_raw_file(points_y_path.c_str(), y, sizeof(G)) != npoints ||
        read_raw_file(points_z_path.c_str(), z, sizeof(G)) != npoints) {
        fprintf(stderr, "Error: Failed to read point coordinates\n");
        delete[] sx;
        delete[] sy;
        delete[] sz;
        delete[] s0;
        delete[] s1;
        delete[] s2;
        delete[] x;
        delete[] y;
        delete[] z;
        return 1;
    }

    printf("Loaded %td surface points, %td surface elements, %td query points\n", 
           nspoints, nselements, npoints);

    // Allocate output array
    T* out = new T[npoints];


    // Initialize output with large distances
    for (ptrdiff_t i = 0; i < npoints; i++) {
      out[i] = std::numeric_limits<T>::max();
    }

    // Compute SDF
    int result = ssdf::sdf(npoints, x, y, z, nselements, s0, s1, s2, nspoints, sx, sy, sz, out);

    if (result != 0) {
        fprintf(stderr, "Error: sdf function returned non-zero: %d\n", result);
        delete[] sx;
        delete[] sy;
        delete[] sz;
        delete[] s0;
        delete[] s1;
        delete[] s2;
        delete[] x;
        delete[] y;
        delete[] z;
        delete[] out;
        return 1;
    }

    // Write output
    FILE* fout = fopen(output_file, "wb");
    if (!fout) {
        fprintf(stderr, "Error: Cannot open output file '%s'\n", output_file);
        delete[] sx;
        delete[] sy;
        delete[] sz;
        delete[] s0;
        delete[] s1;
        delete[] s2;
        delete[] x;
        delete[] y;
        delete[] z;
        delete[] out;
        return 1;
    }

    // Write SDF results (values only, no count)
    if (fwrite(out, sizeof(T), npoints, fout) != npoints) {
        fprintf(stderr, "Error: Failed to write SDF results\n");
        fclose(fout);
        delete[] sx;
        delete[] sy;
        delete[] sz;
        delete[] s0;
        delete[] s1;
        delete[] s2;
        delete[] x;
        delete[] y;
        delete[] z;
        delete[] out;
        return 1;
    }

    fclose(fout);

    // Cleanup
    delete[] sx;
    delete[] sy;
    delete[] sz;
    delete[] s0;
    delete[] s1;
    delete[] s2;
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] out;

    printf("Successfully computed SDF for %td points, written to '%s'\n", npoints, output_file);
    return 0;
}
