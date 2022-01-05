#pragma once

#include "MIComplex.h"

#include <array>
#include <vector>

namespace simplicial_arrangement {

std::array<size_t, 3> mi_cut_3_face(MIComplex<3>& mi_complex,
    size_t cid,
    size_t material_index,
    const std::vector<std::array<size_t, 3>>& subfaces);
} // namespace simplicial_arrangement
