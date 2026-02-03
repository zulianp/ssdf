#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

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
    if (!f) return -1;
    long file_size = get_file_size(f);
    fclose(f);

    if (file_size < 0 || file_size % static_cast<long>(element_size) != 0) return -1;
    return file_size;
}

// Helper function to read raw file and return element count
static ptrdiff_t read_raw_file(const char *filepath, void *data, size_t element_size) {
    FILE *f = fopen(filepath, "rb");
    if (!f) {
        std::fprintf(stderr, "Error: Cannot open file '%s'\n", filepath);
        return -1;
    }

    long file_size = get_file_size(f);
    if (file_size < 0) {
        std::fprintf(stderr, "Error: Cannot determine size of file '%s'\n", filepath);
        std::fclose(f);
        return -1;
    }

    if (file_size % static_cast<long>(element_size) != 0) {
        std::fprintf(stderr,
                     "Error: File '%s' size (%ld) is not a multiple of element size (%zu)\n",
                     filepath,
                     file_size,
                     element_size);
        std::fclose(f);
        return -1;
    }

    ptrdiff_t count = file_size / static_cast<long>(element_size);
    if (std::fread(data, element_size, static_cast<size_t>(count), f) != static_cast<size_t>(count)) {
        std::fprintf(stderr, "Error: Failed to read file '%s'\n", filepath);
        std::fclose(f);
        return -1;
    }

    std::fclose(f);
    return count;
}

template <typename T>
static int write_raw_file(const char *output_file, ptrdiff_t count, const T *data) {
    std::printf("Writing output to '%s'...\n", output_file);
    FILE *fout = std::fopen(output_file, "wb");
    if (!fout) {
        std::fprintf(stderr, "Error: Cannot open output file '%s'\n", output_file);
        return 1;
    }

    if (std::fwrite(data, sizeof(T), static_cast<size_t>(count), fout) != static_cast<size_t>(count)) {
        std::fprintf(stderr, "Error: Failed to write output\n");
        std::fclose(fout);
        return 1;
    }

    std::fclose(fout);
    return 0;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        std::fprintf(stderr,
                     "Usage: %s <udf_file> <inside_u8_file> <output_sdf_file>\n",
                     argv[0]);
        std::fprintf(stderr, "\n");
        std::fprintf(stderr, "Input formats:\n");
        std::fprintf(stderr, "  - udf_file: float udf[npoints] (binary)\n");
        std::fprintf(stderr, "  - inside_u8_file: uint8 inside[npoints] (0=outside, 1=inside)\n");
        std::fprintf(stderr, "\n");
        std::fprintf(stderr, "Output format:\n");
        std::fprintf(stderr, "  - float sdf[npoints] (negative inside)\n");
        return 1;
    }

    const char *udf_file = argv[1];
    const char *inside_file = argv[2];
    const char *output_file = argv[3];

    using T = float;

    const long udf_bytes = check_file_size(udf_file, sizeof(T));
    if (udf_bytes < 0) {
        std::fprintf(stderr, "Error: Cannot open or invalid file '%s'\n", udf_file);
        return 1;
    }

    const long inside_bytes = check_file_size(inside_file, sizeof(uint8_t));
    if (inside_bytes < 0) {
        std::fprintf(stderr, "Error: Cannot open or invalid file '%s'\n", inside_file);
        return 1;
    }

    const ptrdiff_t n_udf = udf_bytes / static_cast<long>(sizeof(T));
    const ptrdiff_t n_inside = inside_bytes / static_cast<long>(sizeof(uint8_t));
    if (n_udf != n_inside) {
        std::fprintf(stderr,
                     "Error: Size mismatch: udf has %td entries, inside has %td entries\n",
                     n_udf,
                     n_inside);
        return 1;
    }

    const ptrdiff_t npoints = n_udf;
    std::printf("Reading %td values...\n", npoints);

    std::vector<T> sdf(npoints);
    std::vector<uint8_t> inside(npoints);
    if (read_raw_file(udf_file, sdf.data(), sizeof(T)) != npoints ||
        read_raw_file(inside_file, inside.data(), sizeof(uint8_t)) != npoints) {
        return 1;
    }

    // Sign: negative inside.
    for (ptrdiff_t i = 0; i < npoints; ++i) {
        if (inside[i]) sdf[i] = -sdf[i];
    }

    return write_raw_file(output_file, npoints, sdf.data());
}

