#pragma once

#include "common.h"

namespace simplicial_arrangement {

// Forward declarations.
template <typename Scalar, int DIM>
class PlaneRepo;

template <int DIM>
struct ARComplex;

namespace internal {

/**
 * Insert a plane into the existing arrangement complex.
 *
 * @param[in] repo            Plane repository.
 * @param[in,out] ar_complex  Current arrangement complex.
 * @param[in] plane_index     The index of the plane to be inserted.
 *
 * @return The index of an existing plane that is coplanar with the inserted
 *         plane if exists.  Otherwise, return `INVALID`.
 */
size_t add_plane(const PlaneRepo<double, 2>& repo, ARComplex<2>& ar_complex, size_t plane_index);
size_t add_plane(const PlaneRepo<Int, 2>& repo, ARComplex<2>& ar_complex, size_t plane_index);
size_t add_plane(const PlaneRepo<double, 3>& repo, ARComplex<3>& ar_complex, size_t plane_index);
size_t add_plane(const PlaneRepo<Int, 3>& repo, ARComplex<3>& ar_complex, size_t plane_index);

} // namespace internal

} // namespace simplicial_arrangement
