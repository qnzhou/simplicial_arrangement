#pragma once

#include <simplicial_arrangement/BSPNode.h>
#include <simplicial_arrangement/common.h>
#include <simplicial_arrangement/cut.h>

#include <array>
#include <type_traits>

namespace simplicial_arrangement {

template <typename Scalar, int DIM>
class SimplicialArrangement
{
public:
    static_assert(DIM == 2 || DIM == 3, "Only 2D and 3D arrangement are supported");
    static_assert(std::is_same<Scalar, double>::value || std::is_same<Scalar, Int>::value,
        "Only double and 128bit int are supported as Scalar.");

public:
    /**
     * Initialize the arrangement induced by a set of planes.
     */
    void initialize(const std::vector<Plane<Scalar, DIM>>& planes)
    {
        set_planes(planes);
        initialize();
    }

    void initialize() {
        if constexpr (DIM == 2) {
            m_root.cell.edges = {0, 1, 2};
        } else {
            m_root.cell.faces = {{0, {1, 2, 3}}, {1, {0, 3, 2}}, {2, {0, 1, 3}}, {3, {0, 2, 1}}};
        }
        m_root.positive = nullptr;
        m_root.negative = nullptr;

        const size_t num_planes = m_planes.size();
        for (size_t i = DIM; i < num_planes; i++) {
            internal::cut(*this, m_root, i);
        }
    }

    void set_planes(const std::vector<Plane<Scalar, DIM>>& planes) {
        const size_t num_cutting_planes = planes.size();
        m_planes.clear();
        m_planes.reserve(num_cutting_planes + DIM + 1);
        if constexpr (DIM == 2) {
            m_planes.push_back({1, 0, 0});
            m_planes.push_back({0, 1, 0});
            m_planes.push_back({0, 0, 1});
        } else {
            m_planes.push_back({1, 0, 0, 0});
            m_planes.push_back({0, 1, 0, 0});
            m_planes.push_back({0, 0, 1, 0});
            m_planes.push_back({0, 0, 0, 1});
        }
        m_planes.insert(m_planes.end(), planes.begin(), planes.end());
    }

    /**
     * Return the set of planes including the initial DIM planes that form the
     * boundary of the simplex.
     */
    const auto& get_planes() const { return m_planes; }

    const BSPNode<DIM>& get_root() const { return m_root; }

private:
    std::vector<Plane<Scalar, DIM>> m_planes;
    BSPNode<DIM> m_root;
};

} // namespace simplicial_arrangement
