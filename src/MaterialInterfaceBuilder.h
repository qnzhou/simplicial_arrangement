#pragma once

#include "BSPNode.h"
#include "DisjointSets.h"
#include "common.h"

#include <absl/container/flat_hash_map.h>

#include <algorithm>
#include <array>
#include <type_traits>

namespace simplicial_arrangement {

template <typename Scalar, int DIM>
class MaterialInterfaceBuilder
{
public:
    static_assert(DIM == 2 || DIM == 3, "Only 2D and 3D material interface are supported");
    static_assert(std::is_same<Scalar, double>::value || std::is_same<Scalar, Int>::value,
        "Only double and 128bit int are supported as Scalar.");

public:
    MaterialInterfaceBuilder(const std::vector<Plane<Scalar, DIM>>& planes)
    {
        set_planes(planes);
    }

    void set_planes(const std::vector<Plane<Scalar, DIM>>& planes)
    {
        //const size_t num_cutting_planes = m_planes.size();
        //if constexpr (DIM == 2) {
        //    m_simplex_planes[0] = {1, 0, 0};
        //    m_simplex_planes[1] = {0, 1, 0};
        //    m_simplex_planes[2] = {0, 0, 1};
        //    register_vertex({1, 2}, 0);
        //    register_vertex({0, 2}, 1);
        //    register_vertex({0, 1}, 2);
        //    m_vertex_count = 3;
        //} else {
        //    m_simplex_planes[0] = {1, 0, 0, 0};
        //    m_simplex_planes[1] = {0, 1, 0, 0};
        //    m_simplex_planes[2] = {0, 0, 1, 0};
        //    m_simplex_planes[3] = {0, 0, 0, 1};
        //    register_vertex({1, 2, 3}, 0);
        //    register_vertex({0, 2, 3}, 1);
        //    register_vertex({0, 1, 3}, 2);
        //    register_vertex({0, 1, 2}, 3);
        //    m_vertex_count = 4;
        //}
        //m_vertex_index_map.reserve(num_cutting_planes * 4);
        //m_coplanar_planes.init(num_cutting_planes + DIM + 1);
    }

    //void register_vertex(Point<DIM> v, size_t index)
    //{
    //    this->sort(v);
    //    auto itr = m_vertex_index_map.find(v);
    //    if (itr == m_vertex_index_map.end()) {
    //        m_vertex_index_map.insert({std::move(v), index});
    //    } else {
    //        assert(itr->second == index);
    //    }
    //}

    //bool has_vertex(Point<DIM> v) const
    //{
    //    this->sort(v);
    //    return m_vertex_index_map.contains(v);
    //}

    //size_t get_vertex_index(Point<DIM> v) const
    //{
    //    this->sort(v);
    //    auto itr = m_vertex_index_map.find(v);
    //    if (itr == m_vertex_index_map.end())
    //        return INVALID;
    //    else
    //        return itr->second;
    //}

    //size_t get_vertex_count() const { return m_vertex_count; }
    //void bump_vertex_count() { m_vertex_count++; }

    //std::vector<Point<DIM>> extract_vertices() const
    //{
    //    std::vector<Point<DIM>> vertices(m_vertex_count);
    //    for (auto& entry : m_vertex_index_map) {
    //        assert(entry.second < m_vertex_count);
    //        vertices[entry.second] = entry.first;
    //    }
    //    return vertices;
    //}

    Arrangement<DIM> extract_material_interface() {
        // TODO.
        return Arrangement<DIM>();
    }

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
    std::vector<Plane<Scalar, DIM+1>> m_planes;
    utils::DisjointSets m_coplanar_planes;
    absl::flat_hash_map<Point<DIM+1>, size_t> m_vertex_index_map;
    size_t m_vertex_count = 0;
};

} // namespace simplicial_arrangement
