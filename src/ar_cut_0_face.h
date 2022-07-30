#pragma once

#include "ARComplex.h"
#include "PlaneRepo.h"
#include "utils.h"

#include <implicit_predicates/implicit_predicates.h>

#ifdef LOOKUP_TABLE_GEN
#error "This header should not be included in lookup table generation mode."
#endif

namespace simplicial_arrangement {

template <typename Scalar, int DIM>
int8_t ar_cut_0_face(const PlaneRepo<Scalar, DIM>& planes,
    const ARComplex<DIM>& ar_complex,
    size_t vid,
    size_t plane_index)
{
    const auto& p = planes.get_plane(plane_index);
    const auto& v = ar_complex.vertices[vid];
    const auto& p0 = planes.get_plane(v[0]);
    const auto& p1 = planes.get_plane(v[1]);

    if constexpr (DIM == 2) {
#ifdef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
        return utils::signof(
            implicit_predicates::orient2d_nonrobust(p0.data(), p1.data(), p.data()));
#else
        return utils::signof(implicit_predicates::orient2d(p0.data(), p1.data(), p.data()));
#endif
    } else {
        const auto& p2 = planes.get_plane(v[2]);
#ifdef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
        return utils::signof(
            implicit_predicates::orient3d_nonrobust(p0.data(), p1.data(), p2.data(), p.data()));
#else
        return utils::signof(
            implicit_predicates::orient3d(p0.data(), p1.data(), p2.data(), p.data()));
#endif
    }
}

} // namespace simplicial_arrangement
