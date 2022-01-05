#pragma once

#include "common.h"

namespace simplicial_arrangement {

// Forward declarations.
template <typename Scalar, int DIM>
class SimplicialArrangementBuilder;

template <int DIM>
struct ARComplex;

Arrangement<2> extract_arrangement(SimplicialArrangementBuilder<double, 2>& builder);
Arrangement<2> extract_arrangement(SimplicialArrangementBuilder<Int, 2>& builder);

Arrangement<3> extract_arrangement(SimplicialArrangementBuilder<double, 3>& builder);
Arrangement<3> extract_arrangement(SimplicialArrangementBuilder<Int, 3>& builder);

Arrangement<2> extract_arrangement(ARComplex<2>&& ar_complex);
Arrangement<3> extract_arrangement(ARComplex<3>&& ar_complex);

} // namespace simplicial_arrangement

