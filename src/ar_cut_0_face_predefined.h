#pragma once

#include "ARComplex.h"
#include "PlaneRepo.h"
#include "utils.h"

#include <implicit_predicates/implicit_predicates.h>

namespace simplicial_arrangement {

extern std::array<std::array<bool, 3>, 4> vertex_signs;
extern std::array<bool, 6> edge_signs;

template <typename Scalar, int DIM>
int8_t ar_cut_0_face(const PlaneRepo<Scalar, DIM>& planes,
    const ARComplex<DIM>& ar_complex,
    size_t vid,
    size_t plane_index)
{
    if constexpr (DIM == 2) {
        throw std::runtime_error("Only 3D lookup table generation is supported.");
    }

    ROBUST_ASSERT(plane_index == 4 || plane_index == 5);
    const auto& v = ar_complex.vertices[vid];

    if (v[0] < 4 && v[1] < 4 && v[2] < 4) {
        // Case 1: v is a tet vertex.
        for (size_t i = 0; i < 4; i++) {
            if (i != v[0] && i != v[1] && i != v[2]) {
                return vertex_signs[i][plane_index - 4] ? 1 : -1;
            }
        }
        ROBUST_ASSERT(false);
    }

    // Case 2: v is on a tet edge.
    ROBUST_ASSERT(plane_index == 5);
    size_t edge_key = INVALID;
    size_t P1 = INVALID, P2 = INVALID;
    for (size_t i = 0; i < 3; i++) {
        if (v[i] < 4) {
            if (P1 == INVALID) {
                P1 = v[i];
            } else {
                P2 = v[i];
            }
        }
    }
    if (P1 > P2) std::swap(P1, P2);

    if (P1 == 0 && P2 == 1) {
        edge_key = 5;
    } else if (P1 == 0 && P2 == 2) {
        edge_key = 4;
    } else if (P1 == 0 && P2 == 3) {
        edge_key = 3;
    } else if (P1 == 1 && P2 == 2) {
        edge_key = 2;
    } else if (P1 == 1 && P2 == 3) {
        edge_key = 1;
    } else if (P1 == 2 && P2 == 3) {
        edge_key = 0;
    }

    return edge_signs[edge_key] ? 1 : -1;
}

} // namespace simplicial_arrangement
