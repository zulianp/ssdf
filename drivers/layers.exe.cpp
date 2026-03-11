#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>

#include "graph/layers.hpp"

namespace fs = std::filesystem;

static long get_file_size(FILE *f) {
    std::fseek(f, 0, SEEK_END);
    long size = std::ftell(f);
    std::fseek(f, 0, SEEK_SET);
    return size;
}

static long check_file_size(const char *filepath, size_t element_size) {
    FILE *f = std::fopen(filepath, "rb");
    if (!f) return -1;
    long file_size = get_file_size(f);
    std::fclose(f);
    if (file_size < 0 || file_size % static_cast<long>(element_size) != 0) return -1;
    return file_size;
}

static ptrdiff_t read_raw_file(const char *filepath, void *data, size_t element_size) {
    FILE *f = std::fopen(filepath, "rb");
    if (!f) {
        std::fprintf(stderr, "Error: Cannot open file '%s'\n", filepath);
        return -1;
    }

    const long file_size = get_file_size(f);
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

    const ptrdiff_t count = file_size / static_cast<long>(element_size);
    if (std::fread(data, element_size, static_cast<size_t>(count), f) != static_cast<size_t>(count)) {
        std::fprintf(stderr, "Error: Failed to read file '%s'\n", filepath);
        std::fclose(f);
        return -1;
    }

    std::fclose(f);
    return count;
}

static int write_u8_file(const char *output_file, ptrdiff_t count, const uint8_t *data) {
    std::printf("Writing output to '%s'...\n", output_file);
    FILE *fout = std::fopen(output_file, "wb");
    if (!fout) {
        std::fprintf(stderr, "Error: Cannot open output file '%s'\n", output_file);
        return 1;
    }
    if (std::fwrite(data, sizeof(uint8_t), static_cast<size_t>(count), fout) != static_cast<size_t>(count)) {
        std::fprintf(stderr, "Error: Failed to write output\n");
        std::fclose(fout);
        return 1;
    }
    std::fclose(fout);
    return 0;
}

struct IndexedFile {
    int idx = -1;
    fs::path path;
};

static bool parse_i_index(const std::string &filename, int &out_idx) {
    // Expect "i<digits>.raw"
    if (filename.size() < 6) return false; // minimal: i0.raw
    if (filename[0] != 'i') return false;
    if (filename.rfind(".raw") != filename.size() - 4) return false;
    int idx = 0;
    size_t pos = 1;
    if (pos >= filename.size() - 4) return false;
    for (; pos < filename.size() - 4; ++pos) {
        const char c = filename[pos];
        if (c < '0' || c > '9') return false;
        idx = idx * 10 + (c - '0');
        if (idx > 1'000'000) return false;
    }
    out_idx = idx;
    return true;
}

int main(int argc, char **argv) {
    if (argc != 5) {
        std::fprintf(stderr,
                     "Usage: %s <mesh_folder> <inside_label> <input_labels.uint8> <output_labels.uint8>\n",
                     argv[0]);
        std::fprintf(stderr, "\n");
        std::fprintf(stderr, "Mesh folder must contain connectivity files i*.raw (e.g., i0.raw, i1.raw, ...).\n");
        std::fprintf(stderr, "nxe is determined by the number of i*.raw files found.\n");
        return 1;
    }

    SSDF_TIMER(layers);

    const fs::path mesh_folder = argv[1];
    const int inside_label_int = std::atoi(argv[2]);
    if (inside_label_int < 0 || inside_label_int > 255) {
        std::fprintf(stderr, "Error: inside_label must be in [0, 255]\n");
        return 1;
    }
    const uint8_t inside_label = static_cast<uint8_t>(inside_label_int);
    const char *input_labels_file = argv[3];
    const char *output_labels_file = argv[4];

    if (!fs::is_directory(mesh_folder)) {
        std::fprintf(stderr, "Error: mesh_folder is not a directory: '%s'\n", mesh_folder.string().c_str());
        return 1;
    }

    // Discover i*.raw files
    std::vector<IndexedFile> index_files;
    for (const auto &entry : fs::directory_iterator(mesh_folder)) {
        if (!entry.is_regular_file()) continue;
        const std::string fname = entry.path().filename().string();
        int idx = -1;
        if (!parse_i_index(fname, idx)) continue;
        index_files.push_back({idx, entry.path()});
    }

    if (index_files.empty()) {
        std::fprintf(stderr, "Error: no connectivity files found in '%s' matching i*.raw\n", mesh_folder.string().c_str());
        return 1;
    }

    // Sort by numeric suffix; nxe = number of discovered files
    std::sort(index_files.begin(), index_files.end(), [](const IndexedFile &a, const IndexedFile &b) { return a.idx < b.idx; });
    const int nxe = static_cast<int>(index_files.size());

    // Read input labels to determine npoints
    const long labels_bytes = check_file_size(input_labels_file, sizeof(uint8_t));
    if (labels_bytes < 0) {
        std::fprintf(stderr, "Error: Cannot open or invalid file '%s'\n", input_labels_file);
        return 1;
    }
    const ptrdiff_t npoints = labels_bytes / static_cast<long>(sizeof(uint8_t));
    if (npoints <= 0) {
        std::fprintf(stderr, "Error: input_labels file is empty\n");
        return 1;
    }

    std::vector<uint8_t> layers(static_cast<size_t>(npoints));
    if (read_raw_file(input_labels_file, layers.data(), sizeof(uint8_t)) != npoints) return 1;

    // Convert labels -> seeds for ssdf::layers(): layer==0 are seeds, everything else is non-zero.
    for (ptrdiff_t i = 0; i < npoints; ++i) {
        layers[i] = (layers[i] == inside_label) ? uint8_t(0) : uint8_t(1);
    }

    // Read connectivity arrays
    using I = int;
    std::vector<std::vector<I>> elems(static_cast<size_t>(nxe));
    ptrdiff_t nelements = -1;

    for (int k = 0; k < nxe; ++k) {
        const fs::path &p = index_files[k].path;
        const long bytes = check_file_size(p.string().c_str(), sizeof(I));
        if (bytes < 0) {
            std::fprintf(stderr, "Error: Cannot open or invalid file '%s'\n", p.string().c_str());
            return 1;
        }
        const ptrdiff_t count = bytes / static_cast<long>(sizeof(I));
        if (nelements < 0) nelements = count;
        if (count != nelements) {
            std::fprintf(stderr,
                         "Error: element count mismatch: '%s' has %td entries, expected %td\n",
                         p.string().c_str(),
                         count,
                         nelements);
            return 1;
        }

        elems[k].resize(static_cast<size_t>(nelements));
        if (read_raw_file(p.string().c_str(), elems[k].data(), sizeof(I)) != nelements) return 1;
    }

    std::vector<I *> elem_ptrs(static_cast<size_t>(nxe));
    for (int k = 0; k < nxe; ++k) elem_ptrs[k] = elems[k].data();

    std::printf("Running layers on %td points, %td elements, nxe=%d (inside_label=%u)\n",
                npoints,
                nelements,
                nxe,
                static_cast<unsigned>(inside_label));

    int SSDF_MAX_LAYERS = 255;
    SSDF_READ_ENV(SSDF_MAX_LAYERS, std::stoi);

    const int rc = ssdf::layers<I>(nelements, nxe, elem_ptrs.data(), npoints, layers.data(), SSDF_MAX_LAYERS);
    if (rc != 0) {
        std::fprintf(stderr, "Error: ssdf::layers returned non-zero: %d\n", rc);
        return 1;
    }

    return write_u8_file(output_labels_file, npoints, layers.data());
}
