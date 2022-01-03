#pragma once

#include "common.h"

#include <array>
#include <vector>

namespace simplicial_arrangement {

template <int DIM>
using ARVertex = Point<DIM>;

template <int DIM>
struct AREdge;

template <int DIM>
struct ARFace;

template <int DIM>
struct ARCell;

template <>
struct AREdge<2>
{
    std::array<size_t, 2> vertices = {INVALID, INVALID}; ///< ordered.
    size_t supporting_plane = INVALID;
};

template <>
struct AREdge<3>
{
    std::array<size_t, 2> vertices = {INVALID, INVALID}; ///< unordered.
    std::array<size_t, 2> supporting_planes = {INVALID, INVALID};
};

template <>
struct ARFace<2>
{
    std::vector<size_t> edges; ///< ordered.
    std::vector<bool> signs; ///< sign of each implicit fucntion.
};

template <>
struct ARFace<3>
{
    std::vector<size_t> edges; ///< ordered.
    size_t supporting_plane = INVALID;
};

template <>
struct ARCell<3>
{
    std::vector<size_t> faces; ///< unordered.
    std::vector<bool> signs; ///< sign of each implicit function.
};

template <int DIM>
struct ARComplex;

template <>
struct ARComplex<2>
{
    std::vector<ARVertex<2>> vertices;
    std::vector<AREdge<2>> edges;
    std::vector<ARFace<2>> faces;
};

template <>
struct ARComplex<3>
{
    std::vector<ARVertex<3>> vertices;
    std::vector<AREdge<3>> edges;
    std::vector<ARFace<3>> faces;
    std::vector<ARCell<3>> cells;
};

} // namespace simplicial_arrangement
