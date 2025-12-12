#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>
#include "ssdf.hpp"

// Helper function to get file size in bytes
static long get_file_size(FILE *f) {
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    fseek(f, 0, SEEK_SET);
    return size;
}

// Helper function to read raw file and return size
static ptrdiff_t read_raw_file(const char *filepath, void *data, size_t element_size) {
    FILE *f = fopen(filepath, "rb");
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
        fprintf(stderr,
                "Error: File '%s' size (%ld) is not a multiple of element size (%zu)\n",
                filepath,
                file_size,
                element_size);
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

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr,
                "Usage: %s <surf_folder1> [surf_folder2] ... <points_folder> "
                "<output_file>\n",
                argv[0]);
        fprintf(stderr, "\n");
        fprintf(stderr, "Surface folder(s) should contain:\n");
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
        fprintf(stderr, "\n");
        fprintf(stderr, "Environment variables:\n");
        fprintf(stderr,
                "  - INPUT_SDF: Path to input SDF file to initialize from "
                "(optional)\n");
        return 1;
    }

    // Parse arguments: last two are points_folder and output_file
    // All others are surface folders
    if (argc < 4) {
        fprintf(stderr,
                "Error: Need at least one surface folder, points folder, "
                "and output file\n");
        return 1;
    }

    const char *points_folder = argv[argc - 2];
    const char *output_file = argv[argc - 1];

    // Collect all surface folders
    std::vector<const char *> surf_folders;
    for (int i = 1; i < argc - 2; i++) {
        surf_folders.push_back(argv[i]);
    }

    using G = float;
    using T = float;
    using I = int;

    // Helper function to construct file paths
    auto make_path = [](const char *folder, const char *filename) -> std::string {
        std::string path = folder;
        if (path.back() != '/') {
            path += '/';
        }
        path += filename;
        return path;
    };

    // Accumulate surface data from all folders
    std::vector<G> sx_vec, sy_vec, sz_vec;
    std::vector<I> s0_vec, s1_vec, s2_vec;

    ptrdiff_t total_nspoints = 0;
    ptrdiff_t total_nselements = 0;

    for (size_t folder_idx = 0; folder_idx < surf_folders.size(); folder_idx++) {
        const char *surf_folder = surf_folders[folder_idx];
        printf("Loading surface folder %zu: %s\n", folder_idx + 1, surf_folder);

        // Read surface point coordinates
        std::string surf_x_path = make_path(surf_folder, "x.raw");
        std::string surf_y_path = make_path(surf_folder, "y.raw");
        std::string surf_z_path = make_path(surf_folder, "z.raw");

        // Determine size by reading one file
        FILE *test_file = fopen(surf_x_path.c_str(), "rb");
        if (!test_file) {
            fprintf(stderr, "Error: Cannot open surface file '%s'\n", surf_x_path.c_str());
            return 1;
        }
        long file_size = get_file_size(test_file);
        fclose(test_file);

        if (file_size < 0 || file_size % sizeof(G) != 0) {
            fprintf(stderr, "Error: Invalid surface x.raw file size in '%s'\n", surf_folder);
            return 1;
        }

        ptrdiff_t nspoints = file_size / sizeof(G);

        // Resize vectors to accommodate new data
        size_t old_size = sx_vec.size();
        sx_vec.resize(old_size + nspoints);
        sy_vec.resize(old_size + nspoints);
        sz_vec.resize(old_size + nspoints);

        if (read_raw_file(surf_x_path.c_str(), &sx_vec[old_size], sizeof(G)) != nspoints ||
            read_raw_file(surf_y_path.c_str(), &sy_vec[old_size], sizeof(G)) != nspoints ||
            read_raw_file(surf_z_path.c_str(), &sz_vec[old_size], sizeof(G)) != nspoints) {
            fprintf(stderr, "Error: Failed to read surface coordinates from '%s'\n", surf_folder);
            return 1;
        }

        // Read surface element indices
        std::string surf_i0_path = make_path(surf_folder, "i0.raw");
        std::string surf_i1_path = make_path(surf_folder, "i1.raw");
        std::string surf_i2_path = make_path(surf_folder, "i2.raw");

        test_file = fopen(surf_i0_path.c_str(), "rb");
        if (!test_file) {
            fprintf(stderr, "Error: Cannot open surface index file '%s'\n", surf_i0_path.c_str());
            return 1;
        }
        file_size = get_file_size(test_file);
        fclose(test_file);

        if (file_size < 0 || file_size % sizeof(I) != 0) {
            fprintf(stderr, "Error: Invalid surface i0.raw file size in '%s'\n", surf_folder);
            return 1;
        }

        ptrdiff_t nselements = file_size / sizeof(I);

        // Resize index vectors
        size_t old_idx_size = s0_vec.size();
        s0_vec.resize(old_idx_size + nselements);
        s1_vec.resize(old_idx_size + nselements);
        s2_vec.resize(old_idx_size + nselements);

        // Adjust indices to point to the accumulated surface points
        I index_offset = total_nspoints;

        if (read_raw_file(surf_i0_path.c_str(), &s0_vec[old_idx_size], sizeof(I)) != nselements ||
            read_raw_file(surf_i1_path.c_str(), &s1_vec[old_idx_size], sizeof(I)) != nselements ||
            read_raw_file(surf_i2_path.c_str(), &s2_vec[old_idx_size], sizeof(I)) != nselements) {
            fprintf(stderr, "Error: Failed to read surface indices from '%s'\n", surf_folder);
            return 1;
        }

        // Adjust indices to account for accumulated points
        for (ptrdiff_t i = 0; i < nselements; i++) {
            s0_vec[old_idx_size + i] += index_offset;
            s1_vec[old_idx_size + i] += index_offset;
            s2_vec[old_idx_size + i] += index_offset;
        }

        total_nspoints += nspoints;
        total_nselements += nselements;
        printf("  Loaded %td surface points, %td surface elements\n", nspoints, nselements);
    }

    // Convert vectors to arrays for the sdf function
    ptrdiff_t nspoints = total_nspoints;
    ptrdiff_t nselements = total_nselements;

    G *sx = sx_vec.data();
    G *sy = sy_vec.data();
    G *sz = sz_vec.data();
    I *s0 = s0_vec.data();
    I *s1 = s1_vec.data();
    I *s2 = s2_vec.data();

    // Read point coordinates
    std::string points_x_path = make_path(points_folder, "x.raw");
    std::string points_y_path = make_path(points_folder, "y.raw");
    std::string points_z_path = make_path(points_folder, "z.raw");

    FILE *test_file = fopen(points_x_path.c_str(), "rb");
    if (!test_file) {
        fprintf(stderr, "Error: Cannot open points file '%s'\n", points_x_path.c_str());
        return 1;
    }
    long file_size = get_file_size(test_file);
    fclose(test_file);

    if (file_size < 0 || file_size % sizeof(G) != 0) {
        fprintf(stderr, "Error: Invalid points x.raw file size\n");
        return 1;
    }

    ptrdiff_t npoints = file_size / sizeof(G);

    G *x = new G[npoints];
    G *y = new G[npoints];
    G *z = new G[npoints];

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

    printf("Total: %td surface points, %td surface elements, %td query points\n", nspoints, nselements, npoints);

    // Read INPUT_SDF environment variable
    const char *INPUT_SDF = getenv("INPUT_SDF");

    // Allocate output array
    T *out = new T[npoints];

    if (INPUT_SDF) {
        // Read SDF field from file
        FILE *sdf_file = fopen(INPUT_SDF, "rb");
        if (!sdf_file) {
            fprintf(stderr, "Error: Cannot open input SDF file '%s'\n", INPUT_SDF);
            delete[] x;
            delete[] y;
            delete[] z;
            delete[] out;
            return 1;
        }

        long sdf_file_size = get_file_size(sdf_file);
        if (sdf_file_size < 0 || sdf_file_size % sizeof(T) != 0) {
            fprintf(stderr, "Error: Invalid input SDF file size (%ld bytes)\n", sdf_file_size);
            fclose(sdf_file);
            delete[] x;
            delete[] y;
            delete[] z;
            delete[] out;
            return 1;
        }

        ptrdiff_t sdf_npoints = sdf_file_size / sizeof(T);
        if (sdf_npoints != npoints) {
            fprintf(stderr, "Error: Input SDF file has %td points, expected %td\n", sdf_npoints, npoints);
            fclose(sdf_file);
            delete[] x;
            delete[] y;
            delete[] z;
            delete[] out;
            return 1;
        }

        if (fread(out, sizeof(T), npoints, sdf_file) != npoints) {
            fprintf(stderr, "Error: Failed to read SDF values from '%s'\n", INPUT_SDF);
            fclose(sdf_file);
            delete[] x;
            delete[] y;
            delete[] z;
            delete[] out;
            return 1;
        }

        fclose(sdf_file);
        printf("Initialized SDF from file: %s\n", INPUT_SDF);
    } else {
        // Initialize output with large distances
        for (ptrdiff_t i = 0; i < npoints; i++) {
            out[i] = std::numeric_limits<T>::max();
        }
    }

    // Compute SDF
    int result = ssdf::edf(npoints, x, y, z, nselements, s0, s1, s2, nspoints, sx, sy, sz, out);

    if (result != 0) {
        fprintf(stderr, "Error: sdf function returned non-zero: %d\n", result);
        delete[] x;
        delete[] y;
        delete[] z;
        delete[] out;
        return 1;
    }

    // Write output
    FILE *fout = fopen(output_file, "wb");
    if (!fout) {
        fprintf(stderr, "Error: Cannot open output file '%s'\n", output_file);
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
        delete[] x;
        delete[] y;
        delete[] z;
        delete[] out;
        return 1;
    }

    fclose(fout);

    // Cleanup (vectors will clean up automatically)
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] out;

    printf("Successfully computed SDF for %td points, written to '%s'\n", npoints, output_file);
    return 0;
}
