#pragma once

namespace simplicial_arrangement {
// Forward declarations.
template <typename Scalar, int DIM>
class SimplicialArrangement;

template <int DIM>
struct BSPNode;
} // namespace simplicial_arrangement

namespace simplicial_arrangement::internal {

void cut(const SimplicialArrangement<double, 2>& arrangement, BSPNode<2>& root, size_t cut_plane);
void cut(const SimplicialArrangement<Int, 2>& arrangement, BSPNode<2>& root, size_t cut_plane);

void cut(const SimplicialArrangement<double, 3>& arrangement, BSPNode<3>& root, size_t cut_plane);
void cut(const SimplicialArrangement<Int, 3>& arrangement, BSPNode<3>& root, size_t cut_plane);

} // namespace simplicial_arrangement::internal
