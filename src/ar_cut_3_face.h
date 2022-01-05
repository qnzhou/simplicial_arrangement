#pragma once

#include "ARComplex.h"

#include <array>
#include <vector>

namespace simplicial_arrangement {

std::array<size_t, 3> ar_cut_3_face(ARComplex<3>& ar_complex,
    size_t cid,
    size_t plane_index,
    const std::vector<std::array<size_t, 3>>& subfaces);
} // namespace simplicial_arrangement
