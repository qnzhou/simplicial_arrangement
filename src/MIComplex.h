#pragma once

#include "common.h"

#include <array>
#include <vector>

namespace simplicial_arrangement {

template <int DIM>
struct MIEdge;

template <int DIM>
struct MIFace;

template <int DIM>
struct MICell;

template <int DIM>
struct MIComplex;

template <>
struct MIEdge<2>
{
    std::array<size_t, 2> vertices{INVALID, INVALID};
    size_t positive_material_label = INVALID;
    size_t negative_material_label = INVALID;
};

template <>
struct MIEdge<3>
{
    std::array<size_t, 2> vertices{INVALID, INVALID};
    std::array<size_t, 3> supporting_materials{INVALID, INVALID, INVALID};
};

template <>
struct MIFace<2>
{
    std::vector<size_t> edges; // ordered.
    size_t material_label = INVALID;
};

template <>
struct MIFace<3>
{
    std::vector<size_t> edges; // ordered
    size_t positive_material_label = INVALID;
    size_t negative_material_label = INVALID;
};

template <>
struct MICell<3>
{
    std::vector<size_t> faces;
    size_t material_label = INVALID;
};

template <>
struct MIComplex<2>
{
    std::vector<Joint<2>> vertices;
    std::vector<MIEdge<2>> edges;
    std::vector<MIFace<2>> faces;
};

template <>
struct MIComplex<3>
{
    std::vector<Joint<3>> vertices;
    std::vector<MIEdge<3>> edges;
    std::vector<MIFace<3>> faces;
    std::vector<MICell<3>> cells;
};

} // namespace simplicial_arrangement
