#pragma once

#include "ARComplex.h"

#include <algorithm>
#include <array>
#include <vector>

namespace simplicial_arrangement {

template <int DIM>
std::array<size_t, 3> ar_cut_2_face(ARComplex<DIM>& ar_complex,
    size_t fid,
    size_t plane_index,
    const std::vector<int8_t>& orientations,
    const std::vector<std::array<size_t, 3>>& subedges)
{
    auto& edges = ar_complex.edges;
    auto& faces = ar_complex.faces;
    edges.reserve(edges.size() + 1);
    faces.reserve(faces.size() + 2);

    auto& f = faces[fid];
    const size_t num_bd_edges = f.edges.size();
    logger().debug("face {} has {} bd edges", fid, num_bd_edges);

    std::vector<size_t> positive_subedges, negative_subedges;
    positive_subedges.reserve(num_bd_edges);
    negative_subedges.reserve(num_bd_edges);

    AREdge<DIM> cut_edge;
    size_t cut_edge_index = INVALID;
    bool face_is_coplanar = true;
    size_t cut_edge_positive_location = INVALID;
    size_t cut_edge_negative_location = INVALID;

    auto get_end_vertex = [&](size_t local_eid) {
        size_t curr_eid = f.edges[local_eid];
        size_t next_eid = f.edges[(local_eid + 1) % num_bd_edges];
        const auto& e0 = edges[curr_eid];
        const auto& e1 = edges[next_eid];
        if (e0.vertices[0] == e1.vertices[0] || e0.vertices[0] == e1.vertices[1]) {
            return e0.vertices[0];
        } else {
            assert(e0.vertices[1] == e1.vertices[0] || e0.vertices[1] == e1.vertices[1]);
            return e0.vertices[1];
        }
    };

    for (size_t j = 0; j < num_bd_edges; j++) {
        const size_t eid = f.edges[j];

        bool last_positive = false;
        bool last_negative = false;
        size_t intersection_point = subedges[eid][2];
        size_t positive_subedge = subedges[eid][0];
        size_t negative_subedge = subedges[eid][1];
        logger().debug("{}: intersection point {}", eid, intersection_point);

        if (positive_subedge == INVALID && negative_subedge == INVALID) {
            // Edge is coplanar with the plane.
            cut_edge_index = eid;
        } else {
            if (positive_subedge != INVALID) {
                positive_subedges.push_back(positive_subedge);
                size_t end_vid = get_end_vertex(j);
                if (orientations[end_vid] <= 0) {
                    // This edge is the last edge with positive subedge.
                    cut_edge_positive_location = positive_subedges.size();
                    last_positive = true;
                }
            }
            if (negative_subedge != INVALID) {
                negative_subedges.push_back(negative_subedge);
                size_t end_vid = get_end_vertex(j);
                if (orientations[end_vid] >= 0) {
                    // This edge is the last edge with negative subedge.
                    cut_edge_negative_location = negative_subedges.size();
                    last_negative = true;
                }
            }
            face_is_coplanar = false;
        }
        if (intersection_point != INVALID) {
            if (last_positive) {
                cut_edge.vertices[0] = intersection_point;
            } else if (last_negative) {
                cut_edge.vertices[1] = intersection_point;
            }
        }
    }
    logger().debug("num positive subedges: {}", positive_subedges.size());
    logger().debug("num negative subedges: {}", negative_subedges.size());
    logger().debug("cut point: {}, {}", cut_edge.vertices[0], cut_edge.vertices[1]);

    if (face_is_coplanar) {
        logger().debug("Face {} is coplanar with cut plane", fid);
        return {INVALID, INVALID, INVALID};
    }

    if (positive_subedges.empty() || negative_subedges.empty()) {
        // No cut.
        if (positive_subedges.empty()) {
            assert(!negative_subedges.empty());
            if constexpr (DIM == 2) f.signs[plane_index] = false;
            if (cut_edge_index != INVALID) {
                // Cut edge is part of the face on the negative side.
                // Need to swap its end points so all cut edges are consistently
                // oriented.
                auto& e = edges[cut_edge_index];
                std::swap(e.vertices[0], e.vertices[1]);
            }
            return {INVALID, fid, cut_edge_index};
        } else {
            assert(!positive_subedges.empty());
            assert(negative_subedges.empty());
            if constexpr (DIM == 2) f.signs[plane_index] = true;
            return {fid, INVALID, cut_edge_index};
        }
    }

    assert(cut_edge_index == INVALID);
    {
        // Insert cut edge.
        if constexpr (DIM == 2) {
            cut_edge.supporting_plane = plane_index;
        } else {
            cut_edge.supporting_planes = {f.supporting_plane, plane_index};
        }
        cut_edge_index = edges.size();
        edges.push_back(std::move(cut_edge));
        logger().debug("Adding cut edge: {}", cut_edge_index);
    }

    // Create subfaces.
    ARFace<DIM> positive_subface, negative_subface;
    if constexpr (DIM == 2) {
        positive_subface.signs = f.signs;
        negative_subface.signs = f.signs;
        positive_subface.signs[plane_index] = true;
        negative_subface.signs[plane_index] = false;
    } else {
        positive_subface.supporting_plane = f.supporting_plane;
        negative_subface.supporting_plane = f.supporting_plane;
    }

    assert(cut_edge_index != INVALID);
    if (cut_edge_positive_location != positive_subedges.size()) {
        std::rotate(positive_subedges.begin(),
            positive_subedges.begin() + cut_edge_positive_location,
            positive_subedges.end());
    }
    if (cut_edge_negative_location != negative_subedges.size()) {
        std::rotate(negative_subedges.begin(),
            negative_subedges.begin() + cut_edge_negative_location,
            negative_subedges.end());
    }
    positive_subedges.push_back(cut_edge_index);
    positive_subface.edges = std::move(positive_subedges);
    assert(positive_subface.edges.size() > 2);
    negative_subedges.push_back(cut_edge_index);
    negative_subface.edges = std::move(negative_subedges);
    assert(negative_subface.edges.size() > 2);

    faces.push_back(std::move(positive_subface));
    faces.push_back(std::move(negative_subface));
    logger().debug("Adding positive subface: {}", faces.size() - 2);
    logger().debug("Adding negative subface: {}", faces.size() - 1);

    // Update face id on each side of cut edge.
    if constexpr (DIM == 2) {
        auto& e = edges[cut_edge_index];
        e.positive_face = faces.size() - 2;
        e.negative_face = faces.size() - 1;
    }

    return {faces.size() - 2, faces.size() - 1, cut_edge_index};
}

} // namespace simplicial_arrangement
