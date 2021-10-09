#pragma once

#include "common.h"

namespace simplicial_arrangement {
// Forward declarations.
template <typename Scalar, int DIM>
class SimplicialArrangementBuilder;

template <int DIM>
struct BSPNode;
} // namespace simplicial_arrangement

namespace simplicial_arrangement::internal {

void cut(SimplicialArrangementBuilder<double, 2>& builder, BSPNode<2>& root, size_t cut_plane);
void cut(SimplicialArrangementBuilder<Int, 2>& builder, BSPNode<2>& root, size_t cut_plane);

void cut(SimplicialArrangementBuilder<double, 3>& builder, BSPNode<3>& root, size_t cut_plane);
void cut(SimplicialArrangementBuilder<Int, 3>& builder, BSPNode<3>& root, size_t cut_plane);

} // namespace simplicial_arrangement::internal
