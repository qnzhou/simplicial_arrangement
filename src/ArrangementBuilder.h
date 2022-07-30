#pragma once

#include "ARComplex.h"
#include "DisjointSets.h"
#include "PlaneRepo.h"
#include "add_plane.h"
#include "common.h"
#include "initialize_simplicial_ar_complex.h"
#include "lookup_table_forward_declarations.h"

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

        bool need_full_compute = true;
        if (use_lookup_table) {
            const auto* ar = lookup(planes);
            if (ar != nullptr) {
                m_arrangement = *ar;
                need_full_compute = false;
            }
        }

        if (need_full_compute) {
            auto ar_complex = initialize_simplicial_ar_complex<DIM>(num_planes + DIM + 1);
            size_t unique_plane_count = 0;
            for (size_t i = 0; i < num_planes; i++) {
                size_t plane_id = i + DIM + 1;
                size_t coplanar_plane_id = internal::add_plane(m_planes, ar_complex, plane_id);
                if (coplanar_plane_id != INVALID) {
                    m_coplanar_planes.merge(plane_id, coplanar_plane_id);
                } else {
                    unique_plane_count++;
                }
            }
            m_arrangement = extract_arrangement(std::move(ar_complex));
            if (unique_plane_count != num_planes) {
                // Only popularize unqiue plane structure if duplicate planes are
                // detected.
                extract_unique_planes();
            }
        }
    }

    const Arrangement<DIM>& get_arrangement() const { return m_arrangement; }
    Arrangement<DIM>& get_arrangement() { return m_arrangement; }
    Arrangement<DIM>&& export_arrangement() && { return std::move(m_arrangement); }

private:
    const Arrangement<DIM>* lookup(const std::vector<Plane<Scalar, DIM>>& planes) const
    {
        extern std::vector<Arrangement<3>> ar_data;
        extern std::vector<size_t> ar_indices;
        const size_t num_planes = planes.size();

        if constexpr (DIM == 3) {
            if (num_planes == 1) {
                const size_t outer_index = ar_compute_outer_index(planes[0]);
                if (outer_index == INVALID) return nullptr;

                logger().debug("AR lookup outer index: {}", outer_index);
                const size_t start_idx = ar_indices[outer_index];
                assert(ar_indices[outer_index + 1] == start_idx + 1);

                logger().debug("AR lookup data index: {}", start_idx);
                return &ar_data[start_idx];
            } else if (num_planes == 2) {
                const size_t outer_index = ar_compute_outer_index(planes[0], planes[1]);
                if (outer_index == INVALID) return nullptr;

                logger().debug("AR lookup outer index: {}", outer_index);
                const size_t start_idx = ar_indices[outer_index];
                const size_t end_idx = ar_indices[outer_index + 1];

                if (end_idx == start_idx + 1) {
                    logger().debug("AR lookup data index: {}", start_idx);
                    return &ar_data[start_idx];
                } else if (end_idx > start_idx) {
                    const size_t inner_index =
                        ar_compute_inner_index(outer_index, planes[0], planes[1]);
                    if (inner_index == INVALID) return nullptr;
                    assert(inner_index < end_idx - start_idx);
                    logger().debug("AR lookup inner index: {}", inner_index);
                    logger().debug("AR lookup data index: {}", start_idx + inner_index);
                    return &ar_data[start_idx + inner_index];
                }
            }
        }
        return nullptr;
    }

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
