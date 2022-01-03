#pragma once

#include "common.h"

namespace simplicial_arrangement {

// Forward declarations.
template <typename Scalar, int DIM>
class PlaneRepo;

template <int DIM>
struct ARComplex;

namespace internal {

void add_plane(const PlaneRepo<double, 2>& repo, ARComplex<2>& ar_complex, size_t plane_index);
void add_plane(const PlaneRepo<Int, 2>& repo, ARComplex<2>& ar_complex, size_t plane_index);

void add_plane(const PlaneRepo<double, 3>& repo, ARComplex<3>& ar_complex, size_t plane_index);
void add_plane(const PlaneRepo<Int, 3>& repo, ARComplex<2>& ar_complex, size_t plane_index);

} // namespace internal

} // namespace simplicial_arrangement
