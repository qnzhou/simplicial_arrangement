#pragma once

#include "common.h"

#include <memory>
#include <vector>

namespace simplicial_arrangement {

template <int DIM>
struct BSPNode
{
    size_t separating_plane = INVALID;
    Cell<DIM> cell;
    std::unique_ptr<BSPNode<DIM>> negative;
    std::unique_ptr<BSPNode<DIM>> positive;
};

} // namespace simplicial_arrangement
