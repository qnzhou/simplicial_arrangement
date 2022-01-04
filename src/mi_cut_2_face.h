#pragma once

#include "MIComplex.h"
#include "utils.h"

#include <array>
#include <vector>

namespace simplicial_arrangement {

template <int DIM>
std::array<size_t, 3> mi_cut_2_face(MIComplex<DIM>& mi_complex,
    size_t fid,
    size_t material_index,
    const std::vector<int8_t>& orientations,
    const std::vector<std::array<size_t, 3>>& subedges)
{
    auto& edges = mi_complex.edges;
    auto& faces = mi_complex.faces;
    edges.reserve(edges.size() + 1);
    faces.reserve(faces.size() + 2);

    const auto& f = faces[fid];
    const size_t num_bd_edges = f.edges.size();
    logger().debug("face {} has {} bd edges", fid, num_bd_edges);
    assert(utils::edges_are_ordered(edges, f.edges));

    std::vector<size_t> positive_subedges, negative_subedges;
    positive_subedges.reserve(num_bd_edges + 1);
    negative_subedges.reserve(num_bd_edges + 1);

    MIEdge<DIM> cut_edge;
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
            // Edge is collinear with material.
            cut_edge_index = eid;
        } else {
            if (positive_subedge != INVALID) {
                positive_subedges.push_back(positive_subedge);
                size_t end_vid = get_end_vertex(j);
                if (orientations[end_vid] <= 0) {
                    cut_edge_positive_location = positive_subedges.size();
                    last_positive = true;
                }
            }
            if (negative_subedge != INVALID) {
                negative_subedges.push_back(negative_subedge);
                size_t end_vid = get_end_vertex(j);
                if (orientations[end_vid] >= 0) {
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
        logger().debug("Face {} is coplanar with material", fid);
        return {INVALID, INVALID, INVALID};
    }

    if (positive_subedges.empty() || negative_subedges.empty()) {
        // No cut.
        if (positive_subedges.empty()) {
            assert(!negative_subedges.empty());
            if constexpr (DIM == 2) {
                faces[fid].material_label = material_index;
            }
            return {INVALID, fid, cut_edge_index};
        } else {
            assert(!positive_subedges.empty());
            assert(negative_subedges.empty());
            return {fid, INVALID, cut_edge_index};
        }
    }

    if (cut_edge_index == INVALID) {
        // Insert cut edge.
        if constexpr (DIM == 2) {
            cut_edge.positive_material_label = f.material_label;
            cut_edge.negative_material_label = material_index;
        } else {
            cut_edge.supporting_materials = {
                f.positive_material_label, f.negative_material_label, material_index};
        }
        cut_edge_index = edges.size();
        edges.push_back(std::move(cut_edge));
        logger().debug("Adding cut edge: {}", cut_edge_index);
    }

    // Create subfaces.
    MIFace<DIM> positive_subface, negative_subface;
    if constexpr (DIM == 2) {
        positive_subface.material_label = f.material_label;
        negative_subface.material_label = material_index;
    } else {
        positive_subface.positive_material_label = f.positive_material_label;
        positive_subface.negative_material_label = f.negative_material_label;
        negative_subface.positive_material_label = f.positive_material_label;
        negative_subface.negative_material_label = f.negative_material_label;
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
    if constexpr (DIM == 2) {
        // Update material labels for 2D edge.
        for (auto eid : negative_subface.edges) {
            auto& e = edges[eid];
            if (e.positive_material_label == f.material_label) {
                e.positive_material_label = material_index;
            } else {
                assert(e.negative_material_label == f.material_label);
                e.negative_material_label = material_index;
            }
        }
    }

    positive_subedges.push_back(cut_edge_index);
    positive_subface.edges = std::move(positive_subedges);
    assert(positive_subface.edges.size() > 2);
    negative_subedges.push_back(cut_edge_index);
    negative_subface.edges = std::move(negative_subedges);
    assert(negative_subface.edges.size() > 2);

    faces.push_back(std::move(positive_subface));
    faces.push_back(std::move(negative_subface));
    size_t positive_subface_id = faces.size() - 2;
    size_t negative_subface_id = faces.size() - 1;
    logger().debug("Adding positive subface: {}", positive_subface_id);
    logger().debug("Adding negative subface: {}", negative_subface_id);

    assert(utils::edges_are_ordered(edges, faces[positive_subface_id].edges));
    assert(utils::edges_are_ordered(edges, faces[negative_subface_id].edges));

    return {positive_subface_id, negative_subface_id, cut_edge_index};
}

} // namespace simplicial_arrangement
