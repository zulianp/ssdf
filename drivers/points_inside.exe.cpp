#include <cstdio>
#include <cstdlib>
#include <cstring>
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

// Helper function to check if file exists and get its size
static long check_file_size(const char *filepath, size_t element_size) {
    FILE *f = fopen(filepath, "rb");
    if (!f) {
        return -1;
    }
    long file_size = get_file_size(f);
    fclose(f);

    if (file_size < 0 || file_size % element_size != 0) {
        return -1;
    }
    return file_size;
}

// Helper function to read raw file and return element count
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
    if (fread(data, element_size, count, f) != static_cast<size_t>(count)) {
        fprintf(stderr, "Error: Failed to read file '%s'\n", filepath);
        fclose(f);
        return -1;
    }

    fclose(f);
    return count;
}

static int write_u8_file(const char *output_file, ptrdiff_t count, const uint8_t *data) {
    printf("Writing output to '%s'...\n", output_file);
    FILE *fout = fopen(output_file, "wb");
    if (!fout) {
        fprintf(stderr, "Error: Cannot open output file '%s'\n", output_file);
        return 1;
    }

    if (fwrite(data, sizeof(uint8_t), static_cast<size_t>(count), fout) != static_cast<size_t>(count)) {
        fprintf(stderr, "Error: Failed to write output\n");
        fclose(fout);
        return 1;
    }

    fclose(fout);
    return 0;
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr,
                "Usage: %s <surf_folder1> [surf_folder2] ... <points_folder> <output_file>\n",
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
        fprintf(stderr, "  - uint8 inside[npoints] (0=outside, 1=inside)\n");
        return 1;
    }

    SSDF_TIMER(points_inside);

    const char *points_folder = argv[argc - 2];
    const char *output_file = argv[argc - 1];

    std::vector<const char *> surf_folders;
    for (int i = 1; i < argc - 2; i++) {
        surf_folders.push_back(argv[i]);
    }

    using G = float;
    using I = int;

    auto make_path = [](const char *folder, const char *filename) -> std::string {
        std::string path = folder;
        if (!path.empty() && path.back() != '/') {
            path += '/';
        }
        path += filename;
        return path;
    };

    // 1) Check all files exist and determine sizes
    printf("Checking files and determining sizes...\n");

    struct SurfaceInfo {
        const char *folder;
        ptrdiff_t nspoints;
        ptrdiff_t nselements;
    };
    std::vector<SurfaceInfo> surface_infos;

    ptrdiff_t total_nspoints = 0;
    ptrdiff_t total_nselements = 0;

    for (size_t folder_idx = 0; folder_idx < surf_folders.size(); folder_idx++) {
        const char *surf_folder = surf_folders[folder_idx];

        std::string surf_x_path = make_path(surf_folder, "x.raw");
        std::string surf_y_path = make_path(surf_folder, "y.raw");
        std::string surf_z_path = make_path(surf_folder, "z.raw");
        std::string surf_i0_path = make_path(surf_folder, "i0.raw");
        std::string surf_i1_path = make_path(surf_folder, "i1.raw");
        std::string surf_i2_path = make_path(surf_folder, "i2.raw");

        long x_size = check_file_size(surf_x_path.c_str(), sizeof(G));
        if (x_size < 0) {
            fprintf(stderr, "Error: Cannot open or invalid file '%s'\n", surf_x_path.c_str());
            return 1;
        }

        if (check_file_size(surf_y_path.c_str(), sizeof(G)) != x_size ||
            check_file_size(surf_z_path.c_str(), sizeof(G)) != x_size) {
            fprintf(stderr, "Error: Surface coordinate files have mismatched sizes in '%s'\n", surf_folder);
            return 1;
        }

        long i0_size = check_file_size(surf_i0_path.c_str(), sizeof(I));
        if (i0_size < 0) {
            fprintf(stderr, "Error: Cannot open or invalid file '%s'\n", surf_i0_path.c_str());
            return 1;
        }

        if (check_file_size(surf_i1_path.c_str(), sizeof(I)) != i0_size ||
            check_file_size(surf_i2_path.c_str(), sizeof(I)) != i0_size) {
            fprintf(stderr, "Error: Surface index files have mismatched sizes in '%s'\n", surf_folder);
            return 1;
        }

        const ptrdiff_t nspoints = x_size / sizeof(G);
        const ptrdiff_t nselements = i0_size / sizeof(I);

        surface_infos.push_back({surf_folder, nspoints, nselements});
        total_nspoints += nspoints;
        total_nselements += nselements;

        printf("  Surface %zu (%s): %td points, %td elements\n", folder_idx + 1, surf_folder, nspoints, nselements);
    }

    // Points folder
    std::string points_x_path = make_path(points_folder, "x.raw");
    std::string points_y_path = make_path(points_folder, "y.raw");
    std::string points_z_path = make_path(points_folder, "z.raw");

    long points_x_size = check_file_size(points_x_path.c_str(), sizeof(G));
    if (points_x_size < 0) {
        fprintf(stderr, "Error: Cannot open or invalid file '%s'\n", points_x_path.c_str());
        return 1;
    }

    if (check_file_size(points_y_path.c_str(), sizeof(G)) != points_x_size ||
        check_file_size(points_z_path.c_str(), sizeof(G)) != points_x_size) {
        fprintf(stderr, "Error: Point coordinate files have mismatched sizes\n");
        return 1;
    }

    const ptrdiff_t npoints = points_x_size / sizeof(G);
    printf("  Points: %td points\n", npoints);

    // 2) Read points
    printf("Reading points...\n");
    std::vector<G> x(npoints), y(npoints), z(npoints);
    if (read_raw_file(points_x_path.c_str(), x.data(), sizeof(G)) != npoints ||
        read_raw_file(points_y_path.c_str(), y.data(), sizeof(G)) != npoints ||
        read_raw_file(points_z_path.c_str(), z.data(), sizeof(G)) != npoints) {
        fprintf(stderr, "Error: Failed to read point coordinates\n");
        return 1;
    }

    // 3) Read and concatenate surfaces (vertex/index offsets)
    printf("Reading %zu surface mesh(es)...\n", surface_infos.size());
    std::vector<G> sx(total_nspoints), sy(total_nspoints), sz(total_nspoints);
    std::vector<I> s0(total_nselements), s1(total_nselements), s2(total_nselements);

    ptrdiff_t vtx_off = 0;
    ptrdiff_t tri_off = 0;
    for (size_t folder_idx = 0; folder_idx < surface_infos.size(); folder_idx++) {
        const SurfaceInfo &info = surface_infos[folder_idx];
        const char *surf_folder = info.folder;
        const ptrdiff_t nspoints = info.nspoints;
        const ptrdiff_t nselements = info.nselements;

        std::string surf_x_path = make_path(surf_folder, "x.raw");
        std::string surf_y_path = make_path(surf_folder, "y.raw");
        std::string surf_z_path = make_path(surf_folder, "z.raw");
        std::string surf_i0_path = make_path(surf_folder, "i0.raw");
        std::string surf_i1_path = make_path(surf_folder, "i1.raw");
        std::string surf_i2_path = make_path(surf_folder, "i2.raw");

        if (read_raw_file(surf_x_path.c_str(), sx.data() + vtx_off, sizeof(G)) != nspoints ||
            read_raw_file(surf_y_path.c_str(), sy.data() + vtx_off, sizeof(G)) != nspoints ||
            read_raw_file(surf_z_path.c_str(), sz.data() + vtx_off, sizeof(G)) != nspoints) {
            fprintf(stderr, "Error: Failed to read surface coordinates from '%s'\n", surf_folder);
            return 1;
        }

        if (read_raw_file(surf_i0_path.c_str(), s0.data() + tri_off, sizeof(I)) != nselements ||
            read_raw_file(surf_i1_path.c_str(), s1.data() + tri_off, sizeof(I)) != nselements ||
            read_raw_file(surf_i2_path.c_str(), s2.data() + tri_off, sizeof(I)) != nselements) {
            fprintf(stderr, "Error: Failed to read surface indices from '%s'\n", surf_folder);
            return 1;
        }

        // Apply vertex offset to local indices.
#pragma omp parallel for
        for (ptrdiff_t i = 0; i < nselements; i++) {
            s0[tri_off + i] += static_cast<I>(vtx_off);
            s1[tri_off + i] += static_cast<I>(vtx_off);
            s2[tri_off + i] += static_cast<I>(vtx_off);
        }

        vtx_off += nspoints;
        tri_off += nselements;
    }

    // 4) Run inside/outside query
    printf("Running points_inside_bvh for %td points (%td vertices, %td triangles)...\n",
           npoints,
           total_nspoints,
           total_nselements);
    std::vector<uint8_t> inside(npoints, 0);

    const int rc = ssdf::points_inside_bvh<G, I>(npoints,
                                                x.data(),
                                                y.data(),
                                                z.data(),
                                                total_nselements,
                                                s0.data(),
                                                s1.data(),
                                                s2.data(),
                                                total_nspoints,
                                                sx.data(),
                                                sy.data(),
                                                sz.data(),
                                                inside.data());
    if (rc != 0) {
        fprintf(stderr, "Error: points_inside_bvh returned non-zero: %d\n", rc);
        return 1;
    }

    // 5) Write output
    if (write_u8_file(output_file, npoints, inside.data()) != 0) return 1;

    printf("Done. Wrote %td inside/outside flags to '%s'\n", npoints, output_file);
    return 0;
}

