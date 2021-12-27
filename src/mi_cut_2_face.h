#pragma once

#include "MIComplex.h"
#include "MaterialInterfaceBuilder.h"

#include <array>
#include <vector>

namespace simplicial_arrangement {

template <typename Scalar, int DIM>
std::array<size_t, 3> mi_cut_2_face(MaterialInterfaceBuilder<Scalar, DIM>& builder,
    MIComplex<DIM>& mi_complex,
    size_t fid,
    size_t material_index,
    const std::vector<int8_t>& orientations,
    const std::vector<std::array<size_t, 3>>& subedges)
{
    auto& edges = mi_complex.edges;
    auto& faces = mi_complex.faces;

    const auto& f = faces[fid];
    const size_t num_bd_edges = f.edges.size();
    logger().debug("face {} has {} bd edges", fid, num_bd_edges);

    // Find edges that straddle the material boundary.
    size_t first_positive_idx = INVALID;
    size_t first_negative_idx = INVALID;
    size_t cut_edge_index = INVALID;
    MIEdge<DIM> cut_edge;

    auto get_shared_vertex = [&](size_t eid_0, size_t eid_1) -> size_t {
        const auto& e0 = edges[eid_0];
        const auto& e1 = edges[eid_1];
        if (e0.vertices[0] == e1.vertices[0] || e0.vertices[0] == e1.vertices[1]) {
            return e0.vertices[0];
        } else if (e0.vertices[1] == e1.vertices[0] || e0.vertices[1] == e1.vertices[1]) {
            return e1.vertices[0];
        } else {
            logger().error("Edge {} and edge {} do not share a vertex", eid_0, eid_1);
            return INVALID;
        }
    };

    for (size_t j = 0; j < num_bd_edges; j++) {
        const size_t eid = f.edges[j];
        const size_t next_eid = f.edges[(j + 1) % num_bd_edges];

        size_t intersection_point = subedges[eid][2];
        size_t positive_subedge = subedges[eid][0];
        size_t negative_subedge = subedges[eid][1];
        logger().debug("{}: intersection point {}", eid, intersection_point);

        if (positive_subedge == INVALID && negative_subedge == INVALID) {
            // Edge is collinear with material.
            cut_edge_index = eid;
        }
        if (intersection_point == INVALID) continue;

        size_t end_vertex = get_shared_vertex(eid, next_eid);
        auto end_orientation = orientations[end_vertex];

        if (positive_subedge != INVALID && end_orientation > 0) {
            assert(first_positive_idx == INVALID);
            first_positive_idx = j;
            cut_edge.vertices[0] = intersection_point;
        }
        if (negative_subedge != INVALID && end_orientation < 0) {
            assert(first_negative_idx == INVALID);
            first_negative_idx = j;
            cut_edge.vertices[1] = intersection_point;
        }
    }
    logger().debug("first positive edge: {}", first_positive_idx);
    logger().debug("first negative edge: {}", first_negative_idx);
    logger().debug("cut point: {}, {}", cut_edge.vertices[0], cut_edge.vertices[1]);

    // Insert cut edge if necessary.
    if (cut_edge.vertices[0] != INVALID && cut_edge.vertices[1] != INVALID) {
        if constexpr (DIM == 2) {
            cut_edge.positive_material_label = f.material_label;
            cut_edge.negative_material_label = material_index;
        } else {
            cut_edge.supporting_materials = {
                f.positive_material_label, f.negative_material_label, material_index};
        }
        cut_edge_index = edges.size();
        edges.push_back(std::move(cut_edge));
    }

    if (first_positive_idx == INVALID || first_negative_idx == INVALID) {
        // No cut.
        logger().debug("No cut!");
        bool on_positive_side = false;
        for (size_t j = 0; j < num_bd_edges; j++) {
            const size_t eid = f.edges[j];
            const auto& e = edges[eid];
            if (orientations[e.vertices[0]] > 0) {
                on_positive_side = true;
                break;
            } else if (orientations[e.vertices[0]] < 0) {
                on_positive_side = false;
                break;
            }
        }
        if (on_positive_side) {
            return {fid, INVALID, cut_edge_index};
        } else {
            return {INVALID, fid, cut_edge_index};
        }
    }

    // Cross cut.
    MIFace<DIM> positive_subface, negative_subface;
    if constexpr (DIM == 2) {
        positive_subface.material_label = f.material_label;
        negative_subface.material_label = material_index;
    } else {
        positive_subface.positive_material_label = f.positive_material_label;
        negative_subface.negative_material_label = f.negative_material_label;
    }
    positive_subface.edges.reserve(num_bd_edges + 1);
    negative_subface.edges.reserve(num_bd_edges + 1);

    positive_subface.edges.push_back(cut_edge_index);
    for (size_t j = 0; j < num_bd_edges; j++) {
        const size_t eid = f.edges[(j + first_positive_idx) % num_bd_edges];
        size_t positive_subedge = subedges[eid][0];
        if (positive_subedge == INVALID) break;
        positive_subface.edges.push_back(positive_subedge);
    }

    negative_subface.edges.push_back(cut_edge_index);
    for (size_t j = 0; j < num_bd_edges; j++) {
        const size_t eid = f.edges[(j + first_negative_idx) % num_bd_edges];
        size_t negative_subedge = subedges[eid][1];
        if (negative_subedge == INVALID) break;
        negative_subface.edges.push_back(negative_subedge);

        if constexpr (DIM == 2) {
            // Update cell id.
            auto& e = edges[negative_subedge];
            if (e.positive_material_label == f.material_label) {
                e.positive_material_label = material_index;
            } else {
                assert(e.negative_material_label == f.material_label);
                e.negative_material_label = material_index;
            }
        }
    }

    assert(cut_edge_index != INVALID);
    assert(positive_subface.edges.size() > 2);
    assert(negative_subface.edges.size() > 2);

    faces.push_back(std::move(positive_subface));
    faces.push_back(std::move(negative_subface));

    return {faces.size() - 2, faces.size() - 1, cut_edge_index};
}

} // namespace simplicial_arrangement
