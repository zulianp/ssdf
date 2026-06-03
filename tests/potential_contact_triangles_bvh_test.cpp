#include "bvh/bvh.hpp"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>

using G = float;
using T = float;
using I = int;
using F = int;

int main() {
    const ptrdiff_t nspoints = 10;
    const G sx[nspoints] = {
        0.0f,
        1.0f,
        0.0f,
        0.2f,
        0.8f,
        0.2f,
        -0.5f,
        3.0f,
        4.0f,
        3.0f,
    };
    const G sy[nspoints] = {
        0.0f,
        0.0f,
        1.0f,
        0.2f,
        0.2f,
        0.8f,
        1.5f,
        0.0f,
        0.0f,
        1.0f,
    };
    const G sz[nspoints] = {
        0.0f,
        0.0f,
        0.0f,
        0.1f,
        0.1f,
        0.1f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
    };

    const ptrdiff_t nselements = 4;
    const I i0[nselements] = {0, 3, 0, 7};
    const I i1[nselements] = {1, 4, 2, 8};
    const I i2[nselements] = {2, 5, 6, 9};

    ptrdiff_t pc_ptr[nselements + 1] = {-7, -7, -7, -7, -7};
    F *pc_idx = nullptr;
    const T extrusion = 0.3f;

    const int err = ssdf::potential_contact_triangles_bvh<G, T, I, F>(
        nselements, i0, i1, i2, nspoints, sx, sy, sz, extrusion, pc_ptr, &pc_idx);

    assert(err == 0);
    assert(pc_idx != nullptr);

    assert(pc_ptr[0] == 0);
    assert(pc_ptr[1] == 1);
    assert(pc_ptr[2] == 2);
    assert(pc_ptr[3] == 2);
    assert(pc_ptr[4] == 2);

    assert(pc_idx[pc_ptr[0]] == 1);
    assert(pc_idx[pc_ptr[1]] == 0);

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
    }

    std::free(pc_idx);
    return 0;
}
