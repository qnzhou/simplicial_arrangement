#pragma once

#include "ARComplex.h"
#include "DisjointSets.h"
#include "PlaneRepo.h"
#include "add_plane.h"
#include "common.h"
#include "initialize_simplicial_ar_complex.h"

#include <algorithm>
#include <array>
#include <type_traits>

namespace simplicial_arrangement {

template <typename Scalar, int DIM>
class ArrangementBuilder
{
public:
    static_assert(DIM == 2 || DIM == 3, "Only 2D and 3D arrangement are supported");
    static_assert(std::is_same<Scalar, double>::value || std::is_same<Scalar, Int>::value,
        "Only double and 128bit int are supported as Scalar.");

public:
    ArrangementBuilder(const std::vector<Plane<Scalar, DIM>>& planes)
        : m_planes(planes)
    {
        const size_t num_planes = m_planes.get_num_planes();
        m_coplanar_planes.init(num_planes + DIM + 1);
        auto ar_complex = initialize_simplicial_ar_complex<DIM>(num_planes + DIM + 1);
        for (size_t i = 0; i < num_planes; i++) {
            size_t plane_id = i + DIM + 1;
            size_t coplanar_plane_id = internal::add_plane(m_planes, ar_complex, plane_id);
            if (coplanar_plane_id != INVALID) {
                m_coplanar_planes.merge(plane_id, coplanar_plane_id);
            }
        }
        m_arrangement = extract_arrangement(std::move(ar_complex));
        extract_unique_planes();
    }

    const Arrangement<DIM>& get_arrangement() const { return m_arrangement; }
    Arrangement<DIM>& get_arrangement() { return m_arrangement; }
    Arrangement<DIM> export_arrangement() { return std::move(m_arrangement); }

private:
    void extract_unique_planes()
    {
        auto is_plane_consistently_oriented = [&](size_t i1, size_t i2) -> bool {
            const auto& p1 = m_planes.get_plane(i1);
            const auto& p2 = m_planes.get_plane(i2);

            for (size_t i = 0; i < DIM + 1; i++) {
                const Scalar v1 = p1[i];
                const Scalar v2 = p2[i];
                if (v1 == 0 && v2 == 0) continue;
                return (v1 > 0 && v2 > 0) || (v1 < 0 && v2 < 0);
            }
            // plane p1 and p2 are consistently 0 over the cell.
            return true;
        };

        auto& r = m_arrangement;
        auto [coplanar_planes, unique_plane_indices] = m_coplanar_planes.extract_disjoint_sets();
        std::vector<bool> coplanar_plane_orientations(unique_plane_indices.size(), true);
        for (const auto& planes : coplanar_planes) {
            const size_t num_planes = planes.size();
            assert(num_planes > 0);
            coplanar_plane_orientations[planes[0]] = true;
            for (size_t i = 1; i < num_planes; i++) {
                coplanar_plane_orientations[planes[i]] =
                    is_plane_consistently_oriented(planes[0], planes[i]);
            }
        }

        r.unique_plane_indices = std::move(unique_plane_indices);
        r.unique_planes = std::move(coplanar_planes);
        r.unique_plane_orientations = std::move(coplanar_plane_orientations);
    }

private:
    PlaneRepo<Scalar, DIM> m_planes;
    utils::DisjointSets m_coplanar_planes;
    Arrangement<DIM> m_arrangement;
};

} // namespace simplicial_arrangement
