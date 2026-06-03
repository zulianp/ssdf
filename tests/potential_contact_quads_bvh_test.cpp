#include "bvh/bvh.hpp"

#include <algorithm>
#include <cassert>
#include <cstdlib>

using G = float;
using T = float;
using I = int;
using F = int;

static bool has_contact(const ptrdiff_t element,
                        const F contact,
                        const ptrdiff_t *const pc_ptr,
                        const F *const pc_idx) {
    return std::find(pc_idx + pc_ptr[element], pc_idx + pc_ptr[element + 1], contact) != pc_idx + pc_ptr[element + 1];
}

int main() {
    const ptrdiff_t nspoints = 20;
    const G sx[nspoints] = {
        0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 3.0f, 4.0f,
        4.0f, 3.0f, 3.0f, 4.0f, 4.0f, 3.0f, 8.0f, 9.0f, 9.0f, 8.0f,
    };
    const G sy[nspoints] = {
        0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f,
        1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f,
    };
    const G sz[nspoints] = {
        0.0f, 0.0f, 0.0f,  0.0f,  0.2f,  0.2f,  0.2f, 0.2f, 0.0f, 0.0f,
        0.5f, 0.0f, 0.75f, 0.75f, 1.25f, 0.75f, 0.0f, 0.0f, 0.0f, 0.0f,
    };

    const ptrdiff_t nselements = 5;
    const I i0[nselements] = {0, 4, 8, 12, 16};
    const I i1[nselements] = {1, 5, 9, 13, 17};
    const I i2[nselements] = {2, 6, 10, 14, 18};
    const I i3[nselements] = {3, 7, 11, 15, 19};

    ptrdiff_t pc_ptr[nselements + 1] = {-7, -7, -7, -7, -7, -7};
    F *pc_idx = nullptr;
    const T extrusion = 0.3f;

    const int err = ssdf::potential_contact_quads_bvh<G, T, I, F>(
        nselements, i0, i1, i2, i3, nspoints, sx, sy, sz, extrusion, pc_ptr, &pc_idx);

    assert(err == 0);
    assert(pc_idx != nullptr);

    assert(pc_ptr[0] == 0);
    assert(pc_ptr[1] == 1);
    assert(pc_ptr[2] == 2);
    assert(pc_ptr[3] == 3);
    assert(pc_ptr[4] == 4);
    assert(pc_ptr[5] == 4);

    assert(has_contact(0, 1, pc_ptr, pc_idx));
    assert(has_contact(1, 0, pc_ptr, pc_idx));
    assert(has_contact(2, 3, pc_ptr, pc_idx));
    assert(has_contact(3, 2, pc_ptr, pc_idx));

    if (false) {
        auto write_binary = [](const char *path, const auto *data, const ptrdiff_t count) {
            FILE *file = std::fopen(path, "wb");
            assert(file != nullptr);
            const size_t written = std::fwrite(data, sizeof(*data), static_cast<size_t>(count), file);
            assert(written == static_cast<size_t>(count));
            assert(std::fclose(file) == 0);
        };

        write_binary("x.float32", sx, nspoints);
        write_binary("y.float32", sy, nspoints);
        write_binary("z.float32", sz, nspoints);
        write_binary("i0.int32", i0, nselements);
        write_binary("i1.int32", i1, nselements);
        write_binary("i2.int32", i2, nselements);
        write_binary("i3.int32", i3, nselements);
    }

    std::free(pc_idx);
    return 0;
}
