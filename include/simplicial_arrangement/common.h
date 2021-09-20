#pragma once

#include <absl/numeric/int128.h>
#include <spdlog/spdlog.h>

#include <array>
#include <limits>
#include <vector>

namespace simplicial_arrangement {

using Int = absl::int128;

constexpr size_t INVALID = std::numeric_limits<size_t>::max();

spdlog::logger& logger();

/**
 * A plane is defined by the barycentric plane equation:
 *     f0 * b0 + f1 * b1 + f2 * b2 + f3 * b3 = 0 // For 3D
 *     f0 * b0 + f1 * b1 + f2 * b2 = 0           // For 2D
 * where the b's are the barycentric variables, and f's are the
 * plane equation coefficients.  We store the f's for each plane.
 */
template <typename Scalar, int DIM>
using Plane = std::array<Scalar, DIM + 1>;

/**
 * A point is represented as the intersection of planes.  We store the index
 * of the plane here.
 */
template <int DIM>
using Point = std::array<size_t, DIM>;

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
