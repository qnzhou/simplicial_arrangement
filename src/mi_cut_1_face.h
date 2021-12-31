#pragma once

#include "MIComplex.h"

#include <implicit_predicates/implicit_predicates.h>

#include <array>
#include <vector>

namespace simplicial_arrangement {

template <int DIM>
std::array<size_t, 3> mi_cut_1_face(
    MIComplex<DIM>& mi_complex,
    size_t eid,
    size_t material_index,
    const std::vector<int8_t>& orientations)
{
    auto& vertices = mi_complex.vertices;
    auto& edges = mi_complex.edges;

    size_t positive_subedge_id = INVALID;
    size_t negative_subedge_id = INVALID;
    size_t intersection_id = INVALID;

    const auto& e = edges[eid];
    const auto& end_points = e.vertices;
    const auto o0 = orientations[end_points[0]];
    const auto o1 = orientations[end_points[1]];
    logger().debug("orientations {} {}", o0, o1);

    if (o0 == 0) intersection_id = end_points[0];
    if (o1 == 0) intersection_id = end_points[1];

    auto compute_intersection_id = [&]() {
        if constexpr (DIM == 2) {
            size_t m0 = e.positive_material_label;
            size_t m1 = e.negative_material_label;
            vertices.push_back({m0, m1, material_index});
            return vertices.size() - 1;
        } else {
            size_t m0 = e.supporting_materials[0];
            size_t m1 = e.supporting_materials[1];
            size_t m2 = e.supporting_materials[2];
            vertices.push_back({m0, m1, m2, material_index});
            return vertices.size() - 1;
        }
    };

    if (o0 == 0 && o1 == 0) {
        // Collinear!
        logger().debug("Edge {} is collinear with material {}", eid, material_index);
    } else if (o0 >= 0 && o1 >= 0) {
        positive_subedge_id = eid;
    } else if (o0 <= 0 && o1 <= 0) {
        negative_subedge_id = eid;
    } else {
        assert(intersection_id == INVALID);
        intersection_id = compute_intersection_id();
        MIEdge<DIM> positive_subedge, negative_subedge;

        if (o0 > 0 && o1 < 0) {
            positive_subedge.vertices = {end_points[0], intersection_id};
            negative_subedge.vertices = {intersection_id, end_points[1]};
        } else {
            assert(o0 < 0);
            assert(o1 > 0);
            negative_subedge.vertices = {end_points[0], intersection_id};
            positive_subedge.vertices = {intersection_id, end_points[1]};
        }

        // Update material labels.
        if constexpr (DIM == 2) {
            positive_subedge.positive_material_label = e.positive_material_label;
            positive_subedge.negative_material_label = e.negative_material_label;
            negative_subedge.positive_material_label = e.positive_material_label;
            negative_subedge.negative_material_label = e.negative_material_label;
        } else {
            positive_subedge.supporting_materials = e.supporting_materials;
            negative_subedge.supporting_materials = e.supporting_materials;
        }

        edges.push_back(std::move(positive_subedge));
        edges.push_back(std::move(negative_subedge));
        positive_subedge_id = edges.size() - 2;
        negative_subedge_id = edges.size() - 1;
    }
    logger().debug("Cut edge {}: {} {} with intersection {}",
        eid,
        positive_subedge_id,
        negative_subedge_id,
        intersection_id);
    return {positive_subedge_id, negative_subedge_id, intersection_id};
}

} // namespace simplicial_arrangement
