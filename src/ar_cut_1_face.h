#pragma once

#include "ARComplex.h"
#include "robust_assert.h"

#include <array>
#include <vector>

namespace simplicial_arrangement {

template <int DIM>
std::array<size_t, 3> ar_cut_1_face(ARComplex<DIM>& ar_complex,
    size_t eid,
    size_t plane_index,
    const std::vector<int8_t>& orientations)
{
    auto& vertices = ar_complex.vertices;
    auto& edges = ar_complex.edges;
    vertices.reserve(vertices.size() + 1);
    edges.reserve(edges.size() + 2);

    size_t positive_subedge_id = INVALID;
    size_t negative_subedge_id = INVALID;
    size_t intersection_id = INVALID;

    const auto& e = edges[eid];
    const auto& end_points = e.vertices;
    const auto o0 = orientations[end_points[0]];
    const auto o1 = orientations[end_points[1]];
    logger().debug("edge {} has end points ({}, {}) with orientations {} {}",
        eid,
        end_points[0],
        end_points[1],
        o0,
        o1);

    if (o0 == 0) intersection_id = end_points[0];
    if (o1 == 0) intersection_id = end_points[1];

    auto compute_intersection_id = [&]() {
        if constexpr (DIM == 2) {
            size_t p0 = e.supporting_plane;
            logger().debug("Adding cut vertex: {}, {}", p0, plane_index);
            vertices.push_back({p0, plane_index});
            return vertices.size() - 1;
        } else {
            size_t p0 = e.supporting_planes[0];
            size_t p1 = e.supporting_planes[1];
            logger().debug("Adding cut vertex: {}, {}, {}", p0, p1, plane_index);
            vertices.push_back({p0, p1, plane_index});
            return vertices.size() - 1;
        }
    };

    if (o0 == 0 && o1 == 0) {
        // Collinear!
        logger().debug("Edge {} is coplanar with plane {}", eid, plane_index);
    } else if (o0 >= 0 && o1 >= 0) {
        positive_subedge_id = eid;
    } else if (o0 <= 0 && o1 <= 0) {
        negative_subedge_id = eid;
    } else {
        ROBUST_ASSERT(intersection_id == INVALID);
        intersection_id = compute_intersection_id();
        AREdge<DIM> positive_subedge, negative_subedge;

        if (o0 > 0 && o1 < 0) {
            positive_subedge.vertices = {end_points[0], intersection_id};
            negative_subedge.vertices = {intersection_id, end_points[1]};
        } else {
            ROBUST_ASSERT(o0 < 0);
            ROBUST_ASSERT(o1 > 0);
            negative_subedge.vertices = {end_points[0], intersection_id};
            positive_subedge.vertices = {intersection_id, end_points[1]};
        }

        // Update supporting materials.
        if constexpr (DIM == 2) {
            positive_subedge.supporting_plane = e.supporting_plane;
            negative_subedge.supporting_plane = e.supporting_plane;
            positive_subedge.positive_face = e.positive_face;
            positive_subedge.negative_face = e.negative_face;
            negative_subedge.positive_face = e.positive_face;
            negative_subedge.negative_face = e.negative_face;
        } else {
            positive_subedge.supporting_planes = e.supporting_planes;
            negative_subedge.supporting_planes = e.supporting_planes;
        }

        edges.push_back(std::move(positive_subedge));
        edges.push_back(std::move(negative_subedge));
        positive_subedge_id = edges.size() - 2;
        negative_subedge_id = edges.size() - 1;
        logger().debug("Adding positive subedge: {}", positive_subedge_id);
        logger().debug("Adding negative subedge: {}", negative_subedge_id);
    }
    logger().debug("Cut edge {}: {} {} with intersection {}",
        eid,
        positive_subedge_id,
        negative_subedge_id,
        intersection_id);
    return {positive_subedge_id, negative_subedge_id, intersection_id};
}

} // namespace simplicial_arrangement
