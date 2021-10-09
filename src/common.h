#pragma once

#include <simplicial_arrangement/simplicial_arrangement.h>
#include <spdlog/spdlog.h>

#include <array>
#include <limits>
#include <vector>

namespace simplicial_arrangement {

constexpr size_t INVALID = std::numeric_limits<size_t>::max();

spdlog::logger& logger();

template <int DIM>
struct Edge;

template <int DIM>
struct Face;

template <int DIM>
struct Cell;


/**
 * A 2D edge is defined by 3 implicit hyperplanes.  The first end point is
 * the intersection of `prev_plane` and `curr_plane`.  The second end point is
 * the intersection of `next_plane` and `curr_plane`.
 */
template <>
struct Edge<2>
{
    Edge(size_t p, size_t c, size_t n)
        : prev_plane(p)
        , curr_plane(c)
        , next_plane(n)
    {}

    size_t prev_plane;
    size_t curr_plane;
    size_t next_plane;
};

template <>
struct Edge<3>
{
    Edge(size_t s, size_t p, size_t c, size_t n)
        : supporting_plane(s)
        , prev_plane(p)
        , curr_plane(c)
        , next_plane(n)
    {}

    size_t supporting_plane;
    size_t prev_plane;
    size_t curr_plane;
    size_t next_plane;
};

template <>
struct Face<3>
{
    size_t supporting_plane;
    std::vector<size_t> edge_planes;
};

template <>
struct Cell<2>
{
    std::vector<size_t> edges;
};

template <>
struct Cell<3>
{
    std::vector<Face<3>> faces;
};

} // namespace simplicial_arrangement
