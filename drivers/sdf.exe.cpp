#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>
#include "ssdf.hpp"

template <typename T>
int write_raw_file(const char *output_file, size_t count, T *data) {
    // Write output
    printf("Writing output to '%s'...\n", output_file);
    FILE *fout = fopen(output_file, "wb");
    if (!fout) {
        fprintf(stderr, "Error: Cannot open output file '%s'\n", output_file);
        return 1;
    }

    // Write SDF results (values only, no count)
    if (fwrite(data, sizeof(T), count, fout) != count) {
        fprintf(stderr, "Error: Failed to write SDF results\n");
        fclose(fout);
        return 1;
    }

    fclose(fout);
    return 0;
}

template <typename I, typename G>
void export_intervals(const ptrdiff_t nselements,
                      const I *const SSDF_RESTRICT s0,
                      const I *const SSDF_RESTRICT s1,
                      const I *const SSDF_RESTRICT s2,
                      const ptrdiff_t nspoints,
                      const G *const SSDF_RESTRICT sx,
                      const G *const SSDF_RESTRICT sy,
                      const G *const SSDF_RESTRICT sz) {
    // Triangle AABBs
    std::vector<G> ix(nselements), iy(nselements), iz(nselements);

#pragma omp parallel
    {
#pragma omp for
        for (ptrdiff_t i = 0; i < nselements; ++i) {
            const I i0 = s0[i], i1 = s1[i], i2 = s2[i];
            const G x0 = sx[i0], x1 = sx[i1], x2 = sx[i2];
            const G y0 = sy[i0], y1 = sy[i1], y2 = sy[i2];
            const G z0 = sz[i0], z1 = sz[i1], z2 = sz[i2];

            const G xmin = std::min({x0, x1, x2});
            const G xmax = std::max({x0, x1, x2});
            const G ymin = std::min({y0, y1, y2});
            const G ymax = std::max({y0, y1, y2});
            const G zmin = std::min({z0, z1, z2});
            const G zmax = std::max({z0, z1, z2});

            ix[i] = xmax - xmin;
            iy[i] = ymax - ymin;
            iz[i] = zmax - zmin;
        }
    }

    write_raw_file("ix.f32", nselements, ix.data());
    write_raw_file("iy.f32", nselements, iy.data());
    write_raw_file("iz.f32", nselements, iz.data());
}

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
                "  - SSDF_INPUT_SDF: Path to input SDF file to initialize from "
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

    // 1) Check that all files exist and determine sizes
    printf("Checking files and determining sizes...\n");

    struct SurfaceInfo {
        const char *folder;
        ptrdiff_t nspoints;
        ptrdiff_t nselements;
    };
    std::vector<SurfaceInfo> surface_infos;

    ptrdiff_t max_nspoints = 0;
    ptrdiff_t max_nselements = 0;

    for (size_t folder_idx = 0; folder_idx < surf_folders.size(); folder_idx++) {
        const char *surf_folder = surf_folders[folder_idx];

        std::string surf_x_path = make_path(surf_folder, "x.raw");
        std::string surf_y_path = make_path(surf_folder, "y.raw");
        std::string surf_z_path = make_path(surf_folder, "z.raw");
        std::string surf_i0_path = make_path(surf_folder, "i0.raw");
        std::string surf_i1_path = make_path(surf_folder, "i1.raw");
        std::string surf_i2_path = make_path(surf_folder, "i2.raw");

        // Check all files exist
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

        ptrdiff_t nspoints = x_size / sizeof(G);
        ptrdiff_t nselements = i0_size / sizeof(I);

        surface_infos.push_back({surf_folder, nspoints, nselements});

        max_nspoints = std::max(max_nspoints, nspoints);
        max_nselements = std::max(max_nselements, nselements);

        printf("  Surface %zu (%s): %td points, %td elements\n", folder_idx + 1, surf_folder, nspoints, nselements);
    }

    // Check points folder files
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

    ptrdiff_t npoints = points_x_size / sizeof(G);
    printf("  Points: %td points\n", npoints);

    // 2) Read points folder
    printf("Reading points...\n");
    std::vector<G> x(npoints);
    std::vector<G> y(npoints);
    std::vector<G> z(npoints);

    if (read_raw_file(points_x_path.c_str(), x.data(), sizeof(G)) != npoints ||
        read_raw_file(points_y_path.c_str(), y.data(), sizeof(G)) != npoints ||
        read_raw_file(points_z_path.c_str(), z.data(), sizeof(G)) != npoints) {
        fprintf(stderr, "Error: Failed to read point coordinates\n");
        return 1;
    }

    // 3) Read optional input SDF
    const char *SSDF_INPUT_SDF = getenv("SSDF_INPUT_SDF");
    std::vector<T> out(npoints);

    if (SSDF_INPUT_SDF) {
        printf("Reading input SDF from '%s'...\n", SSDF_INPUT_SDF);
        FILE *sdf_file = fopen(SSDF_INPUT_SDF, "rb");
        if (!sdf_file) {
            fprintf(stderr, "Error: Cannot open input SDF file '%s'\n", SSDF_INPUT_SDF);
            return 1;
        }

        long sdf_file_size = get_file_size(sdf_file);
        if (sdf_file_size < 0 || sdf_file_size % sizeof(T) != 0) {
            fprintf(stderr, "Error: Invalid input SDF file size (%ld bytes)\n", sdf_file_size);
            fclose(sdf_file);
            return 1;
        }

        ptrdiff_t sdf_npoints = sdf_file_size / sizeof(T);
        if (sdf_npoints != npoints) {
            fprintf(stderr, "Error: Input SDF file has %td points, expected %td\n", sdf_npoints, npoints);
            fclose(sdf_file);
            return 1;
        }

        if (fread(out.data(), sizeof(T), npoints, sdf_file) != npoints) {
            fprintf(stderr, "Error: Failed to read SDF values from '%s'\n", SSDF_INPUT_SDF);
            fclose(sdf_file);
            return 1;
        }

        fclose(sdf_file);
    } else {
        // Initialize with points box (case where surface points coincide with query points)
        int SSDF_INIT_WITH_VOLUME_BOX = 0;
        SSDF_READ_ENV(SSDF_INIT_WITH_VOLUME_BOX, std::stoi);
        if (SSDF_INIT_WITH_VOLUME_BOX) {
            G xmin, xmax, ymin, ymax, zmin, zmax;
            ssdf::compute_aabb(npoints, x.data(), y.data(), z.data(), &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);
            ssdf::all_points_aabb_signed_distance<G, T>(
                npoints, x.data(), y.data(), z.data(), xmin, xmax, ymin, ymax, zmin, zmax, out.data());
#pragma omp parallel for
            for (ptrdiff_t i = 0; i < npoints; i++) {
                out[i] = std::abs(out[i]);
            }
        } else {
            // Initialize output with large distances
#pragma omp parallel for
            for (ptrdiff_t i = 0; i < npoints; i++) {
                out[i] = std::numeric_limits<T>::max();
            }
        }
    }

    // 4) Allocate reusable surface memory (use maximum size needed)
    printf("Allocating surface memory (max: %td points, %td elements)...\n", max_nspoints, max_nselements);
    std::vector<G> sx(max_nspoints);
    std::vector<G> sy(max_nspoints);
    std::vector<G> sz(max_nspoints);
    std::vector<I> s0(max_nselements);
    std::vector<I> s1(max_nselements);
    std::vector<I> s2(max_nselements);

    // 5) Loop through surfaces, read one at a time and call edf
    printf("Processing surfaces...\n");
    for (size_t folder_idx = 0; folder_idx < surface_infos.size(); folder_idx++) {
        const SurfaceInfo &info = surface_infos[folder_idx];
        const char *surf_folder = info.folder;
        ptrdiff_t nspoints = info.nspoints;
        ptrdiff_t nselements = info.nselements;

        printf("  Processing surface %zu/%zu: %s (%td points, %td elements)\n",
               folder_idx + 1,
               surface_infos.size(),
               surf_folder,
               nspoints,
               nselements);

        // Read surface point coordinates
        std::string surf_x_path = make_path(surf_folder, "x.raw");
        std::string surf_y_path = make_path(surf_folder, "y.raw");
        std::string surf_z_path = make_path(surf_folder, "z.raw");

        if (read_raw_file(surf_x_path.c_str(), sx.data(), sizeof(G)) != nspoints ||
            read_raw_file(surf_y_path.c_str(), sy.data(), sizeof(G)) != nspoints ||
            read_raw_file(surf_z_path.c_str(), sz.data(), sizeof(G)) != nspoints) {
            fprintf(stderr, "Error: Failed to read surface coordinates from '%s'\n", surf_folder);
            return 1;
        }

        // Read surface element indices
        std::string surf_i0_path = make_path(surf_folder, "i0.raw");
        std::string surf_i1_path = make_path(surf_folder, "i1.raw");
        std::string surf_i2_path = make_path(surf_folder, "i2.raw");

        if (read_raw_file(surf_i0_path.c_str(), s0.data(), sizeof(I)) != nselements ||
            read_raw_file(surf_i1_path.c_str(), s1.data(), sizeof(I)) != nselements ||
            read_raw_file(surf_i2_path.c_str(), s2.data(), sizeof(I)) != nselements) {
            fprintf(stderr, "Error: Failed to read surface indices from '%s'\n", surf_folder);
            return 1;
        }

        // export_intervals(nselements, s0.data(), s1.data(), s2.data(), nspoints, sx.data(), sy.data(), sz.data());

        int SSDF_DOUBLE_PRECISION = 0;
        SSDF_READ_ENV(SSDF_DOUBLE_PRECISION, std::stoi);

        double tick = ssdf::time_ms();

        int result = 0;
        if (SSDF_DOUBLE_PRECISION == 1 && !std::is_same<T, double>::value) {
            // Convert input and out to double
            std::vector<double> x_double(npoints);
            std::vector<double> y_double(npoints);
            std::vector<double> z_double(npoints);
            std::vector<double> out_double(npoints);
            std::vector<double> sx_double(sx.size());
            std::vector<double> sy_double(sy.size());
            std::vector<double> sz_double(sz.size());

            for (ptrdiff_t i = 0; i < npoints; ++i) {
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

            result = ssdf::edf_select<double, double, int>(npoints,
                                                           x_double.data(),
                                                           y_double.data(),
                                                           z_double.data(),
                                                           nselements,
                                                           s0.data(),
                                                           s1.data(),
                                                           s2.data(),
                                                           nspoints,
                                                           sx_double.data(),
                                                           sy_double.data(),
                                                           sz_double.data(),
                                                           out_double.data());

            // Convert output to float
            for (ptrdiff_t i = 0; i < npoints; ++i) {
                out[i] = T(out_double[i]);
            }

        } else {
            // Call edf for this surface
            result = ssdf::edf_select(npoints,
                                      x.data(),
                                      y.data(),
                                      z.data(),
                                      nselements,
                                      s0.data(),
                                      s1.data(),
                                      s2.data(),
                                      nspoints,
                                      sx.data(),
                                      sy.data(),
                                      sz.data(),
                                      out.data());
        }

        double tock = ssdf::time_ms();
        printf("Time taken: %f ms\n", tock - tick);

        if (result != 0) {
            fprintf(stderr, "Error: edf function returned non-zero: %d\n", result);
            return 1;
        }
    }

    // Write output
    printf("Writing output to '%s'...\n", output_file);
    FILE *fout = fopen(output_file, "wb");
    if (!fout) {
        fprintf(stderr, "Error: Cannot open output file '%s'\n", output_file);
        return 1;
    }

    // Write SDF results (values only, no count)
    if (fwrite(out.data(), sizeof(T), npoints, fout) != npoints) {
        fprintf(stderr, "Error: Failed to write SDF results\n");
        fclose(fout);
        return 1;
    }

    fclose(fout);

    printf("Successfully computed SDF for %td points across %zu surfaces, written to '%s'\n",
           npoints,
           surface_infos.size(),
           output_file);
    return 0;
}
