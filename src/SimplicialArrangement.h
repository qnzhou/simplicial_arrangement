#pragma once

#include "BSPNode.h"
#include "DisjointSets.h"
#include "common.h"
#include "cut.h"
#include "extract_arrangement.h"

#include <absl/container/flat_hash_map.h>

#include <algorithm>
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
    SimplicialArrangement(const std::vector<Plane<Scalar, DIM>>& planes) {
        initialize(planes);
    }

    /**
     * Initialize the arrangement induced by a set of planes.
     */
    void initialize(const std::vector<Plane<Scalar, DIM>>& planes)
    {
        set_planes(planes);
        initialize();
    }

    void initialize()
    {
        if constexpr (DIM == 2) {
            m_root.cell.edges = {0, 1, 2};
        } else {
            m_root.cell.faces = {{0, {1, 2, 3}}, {1, {0, 3, 2}}, {2, {0, 1, 3}}, {3, {0, 2, 1}}};
        }
        m_root.positive = nullptr;
        m_root.negative = nullptr;

        const size_t num_planes = m_planes.size();
        for (size_t i = DIM + 1; i < num_planes; i++) {
            internal::cut(*this, m_root, i);
        }
    }

    void set_planes(const std::vector<Plane<Scalar, DIM>>& planes)
    {
        const size_t num_cutting_planes = planes.size();
        m_planes.clear();
        m_planes.reserve(num_cutting_planes + DIM + 1);
        if constexpr (DIM == 2) {
            m_planes.push_back({1, 0, 0});
            m_planes.push_back({0, 1, 0});
            m_planes.push_back({0, 0, 1});
            register_vertex({1, 2}, 0);
            register_vertex({0, 2}, 1);
            register_vertex({0, 1}, 2);
            m_vertex_count = 3;
        } else {
            m_planes.push_back({1, 0, 0, 0});
            m_planes.push_back({0, 1, 0, 0});
            m_planes.push_back({0, 0, 1, 0});
            m_planes.push_back({0, 0, 0, 1});
            register_vertex({1, 2, 3}, 0);
            register_vertex({0, 2, 3}, 1);
            register_vertex({0, 1, 3}, 2);
            register_vertex({0, 1, 2}, 3);
            m_vertex_count = 4;
        }
        m_planes.insert(m_planes.end(), planes.begin(), planes.end());
        m_vertex_index_map.reserve(num_cutting_planes * 4);
        m_coplanar_planes.init(num_cutting_planes + DIM + 1);
    }

    /**
     * Return the set of planes including the initial DIM planes that form the
     * boundary of the simplex.
     */
    const auto& get_planes() const { return m_planes; }

    size_t get_num_planes() const { return m_planes.size(); }

    const BSPNode<DIM>& get_root() const { return m_root; }

    void register_vertex(Point<DIM> v, size_t index)
    {
        this->sort(v);
        auto itr = m_vertex_index_map.find(v);
        if (itr == m_vertex_index_map.end()) {
            m_vertex_index_map.insert({std::move(v), index});
        } else {
            assert(itr->second == index);
        }
    }

    bool has_vertex(Point<DIM> v) const
    {
        this->sort(v);
        return m_vertex_index_map.contains(v);
    }

    size_t get_vertex_index(Point<DIM> v) const
    {
        this->sort(v);
        auto itr = m_vertex_index_map.find(v);
        if (itr == m_vertex_index_map.end())
            return INVALID;
        else
            return itr->second;
    }

    size_t get_vertex_count() const { return m_vertex_count; }
    void bump_vertex_count() { m_vertex_count++; }

    std::vector<Point<DIM>> extract_vertices() const
    {
        std::vector<Point<DIM>> vertices(m_vertex_count);
        for (auto& entry : m_vertex_index_map) {
            assert(entry.second < m_vertex_count);
            vertices[entry.second] = entry.first;
        }
        return vertices;
    }

    void register_coplanar_planes(size_t p0, size_t p1) { m_coplanar_planes.merge(p0, p1); }

    std::tuple<std::vector<std::vector<size_t>>, std::vector<size_t>> extract_coplanar_planes()
    {
        return m_coplanar_planes.extract_disjoint_sets();
    }

    Arrangement<DIM> extract_arrangement() { return internal::extract_arrangement(*this); }

private:
    inline void sort(Point<DIM>& p) const
    {
        if constexpr (DIM == 2) {
            if (p[0] > p[1]) std::swap(p[0], p[1]);
        } else if (DIM == 3) {
            if (p[0] > p[1]) std::swap(p[0], p[1]);
            if (p[0] > p[2]) std::swap(p[0], p[2]);
            if (p[1] > p[2]) std::swap(p[1], p[2]);
        }
    }

private:
    std::vector<Plane<Scalar, DIM>> m_planes;
    utils::DisjointSets m_coplanar_planes;
    absl::flat_hash_map<Point<DIM>, size_t> m_vertex_index_map;
    BSPNode<DIM> m_root;
    size_t m_vertex_count = 0;
};

} // namespace simplicial_arrangement
