#pragma once

#include "common.h"

namespace simplicial_arrangement {
// Forward declarations.
template <typename Scalar, int DIM>
class SimplicialArrangementBuilder;
} // namespace simplicial_arrangement

namespace simplicial_arrangement::internal {

Arrangement<2> extract_arrangement(SimplicialArrangementBuilder<double, 2>& builder);
Arrangement<2> extract_arrangement(SimplicialArrangementBuilder<Int, 2>& builder);

Arrangement<3> extract_arrangement(SimplicialArrangementBuilder<double, 3>& builder);
Arrangement<3> extract_arrangement(SimplicialArrangementBuilder<Int, 3>& builder);

} // namespace simplicial_arrangement::internal

