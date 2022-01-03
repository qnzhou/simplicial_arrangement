#pragma once

#include "common.h"

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
    ArrangementBuilder(const std::vector<Plane<Scalar, DIM>>& planes) : m_planes(planes) {
        initialize();
    }

    void initialize()
    {
        const size_t num_planes = m_planes.get_num_planes() + DIM + 1;
        auto ar_complex = initialize_simplicial_ar_complex<DIM>();
        for (size_t i = 0; i < num_planes; i++) {
            // TODO
            // internal::add_plane(*this, ar_complex, i + DIM + 1);
        }
    }

    const Plane<Scalar, DIM>& get_plane(size_t i) const {
        return m_planes.get_plane(i);
    }

    void register_coplanar_planes(size_t p0, size_t p1) { m_coplanar_planes.merge(p0, p1); }

    std::tuple<std::vector<std::vector<size_t>>, std::vector<size_t>> extract_coplanar_planes()
    {
        return m_coplanar_planes.extract_disjoint_sets();
    }

    const Arrangement<DIM>& get_arrangement() const { return m_arrangement; }
    Arrangement<DIM>& get_arrangement() { return m_arrangement; }

private:
    PlaneRepo<Scalar, DIM> m_planes;
    utils::DisjointSets m_coplanar_planes;
    Arrangement<DIM> m_arrangement;
};

} // namespace simplicial_arrangement
