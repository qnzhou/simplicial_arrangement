#pragma once

#include "MIComplex.h"
#include "MaterialRepo.h"
#include "utils.h"

#include <vector>

namespace simplicial_arrangement {

extern std::array<std::array<bool, 3>, 4> vertex_signs; // { m0 > m1, m1 > m2, m0 > m2 } x 4.
extern std::array<bool, 6> edge_signs; // intersection of m0 and m1 > m2 for each edge.

/**
 * This method is only used for enumerating the lookup table. Do not use in
 * practice.
 */
template <typename Scalar>
int8_t mi_cut_0_face([[maybe_unused]] const MaterialRepo<Scalar, 3>& materials,
    const MIComplex<3>& mi_complex,
    size_t vid,
    size_t material_index)
{
    assert(material_index > 3);
    const auto& vertices = mi_complex.vertices;
    const auto& p = vertices[vid];

    size_t m0 = INVALID, m1 = INVALID;
    for (size_t i = 0; i < 4; i++) {
        if (p[i] <= 3) continue;
        if (m0 == INVALID) {
            m0 = p[i];
        } else if (m1 == INVALID) {
            m1 = p[i];
        } else {
            logger().error(
                "Query point ({}, {}, {}, {}) is not on a vertex or edge!", p[0], p[1], p[2], p[3]);
        }
    }
    assert(m0 != INVALID);
    assert(m0 < material_index);

    // Case 1: p is a tet vertex.
    if (m1 == INVALID) {
        assert(m0 < material_index);
        for (size_t i = 0; i < 4; i++) {
            if (i != p[0] && i != p[1] && i != p[2] && i != p[3]) {
                if (m0 == 4) {
                    if (material_index == 5) {
                        return vertex_signs[i][0] ? 1 : -1;
                    } else {
                        return vertex_signs[i][2] ? 1 : -1;
                    }
                } else {
                    return vertex_signs[i][1] ? 1 : -1;
                }
            }
        }
    }

    // Case 2: p is on a tet edge.
    if (m0 > m1) std::swap(m0, m1);
    assert(m1 < material_index);

    size_t edge_key = INVALID;

    size_t b0 = INVALID, b1 = INVALID;
    for (size_t i = 0; i < 4; i++) {
        if (p[i] <= 3) {
            if (b0 == INVALID) {
                b0 = p[i];
            } else {
                assert(b1 == INVALID);
                b1 = p[i];
            }
        }
    }
    if (b0 > b1) std::swap(b0, b1);

    if (b0 == 0 && b1 == 1) {
        edge_key = 5;
    } else if (b0 == 0 && b1 == 2) {
        edge_key = 4;
    } else if (b0 == 0 && b1 == 3) {
        edge_key = 3;
    } else if (b0 == 1 && b1 == 2) {
        edge_key = 2;
    } else if (b0 == 1 && b1 == 3) {
        edge_key = 1;
    } else if (b0 == 2 && b1 == 3) {
        edge_key = 0;
    }

    logger().debug("mi_cut_0d: {}, {}, {} at edge {} -> {}", m0, m1, material_index, edge_key, edge_signs[edge_key]);
    return edge_signs[edge_key] ? 1 : -1;
}

template <typename Scalar>
int8_t mi_cut_0_face([[maybe_unused]] const MaterialRepo<Scalar, 2>& materials,
    [[maybe_unused]] const MIComplex<2>& mi_complex,
    [[maybe_unused]] size_t vid,
    [[maybe_unused]] size_t material_index)
{
    return 0;
}

} // namespace simplicial_arrangement
