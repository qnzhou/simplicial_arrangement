#pragma once

#include "common.h"

namespace simplicial_arrangement {
// Forward declarations.
template <typename Scalar, int DIM>
class SimplicialArrangement;
} // namespace simplicial_arrangement

namespace simplicial_arrangement::internal {

Arrangement<2> extract_arrangement(SimplicialArrangement<double, 2>& arrangement);
Arrangement<2> extract_arrangement(SimplicialArrangement<Int, 2>& arrangement);

Arrangement<3> extract_arrangement(SimplicialArrangement<double, 3>& arrangement);
Arrangement<3> extract_arrangement(SimplicialArrangement<Int, 3>& arrangement);

} // namespace simplicial_arrangement::internal

