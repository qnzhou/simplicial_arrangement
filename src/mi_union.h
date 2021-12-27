#pragma once

#include "MIComplex.h"
#include "MaterialInterfaceBuilder.h"

#include <array>
#include <vector>

namespace simplicial_arrangement {

template <typename Scalar>
size_t mi_union_negative_faces(MaterialInterfaceBuilder<Scalar, 2>& builder,
    MIComplex<2>& mi_complex,
    size_t material_index,
    const std::vector<std::array<size_t, 3>>& subfaces)
{
    const size_t num_negative_faces = std::count_if(
        subfaces.begin(), subfaces.end(), [](const auto& val) { return val[1] != INVALID; });
    if (num_negative_faces == 0) return INVALID;
    logger().debug("combining {} negative cells.", num_negative_faces);

    auto& vertices = mi_complex.vertices;
    auto& edges = mi_complex.edges;
    auto& faces = mi_complex.faces;

    MIFace<2> combined_face;
    combined_face.material_label = material_index;
    const size_t total_num_vertices = mi_complex.vertices.size();
    const size_t total_num_edges = mi_complex.edges.size();

    std::vector<bool> on_border(total_num_edges, false);
    std::vector<size_t> next_edge(total_num_vertices, INVALID);

    for (const auto& subface : subfaces) {
        const auto fid = subface[1];
        if (fid == INVALID) continue;
        const auto& f = faces[fid];
        for (const auto& eid : f.edges) {
            on_border[eid] = !on_border[eid];
        }
    }

    auto get_ordered_vertex = [&](size_t eid, size_t i) {
        assert(i < 2);
        const auto& e = edges[eid];
        if (e.positive_material_label == material_index) {
            return e.vertices[(i + 1) % 2];
        } else {
            assert(e.negative_material_label == material_index);
            return e.vertices[i];
        }
    };

    size_t num_bd_edges = 0;
    size_t start_edge = INVALID;
    for (size_t i = 0; i < total_num_edges; i++) {
        if (!on_border[i]) continue;
        start_edge = i;
        next_edge[get_ordered_vertex(i, 0)] = i;
        num_bd_edges++;
    }

    combined_face.edges.reserve(num_bd_edges);
    size_t eid = start_edge;
    for (size_t i = 0; i < num_bd_edges && eid != INVALID; i++) {
        combined_face.edges.push_back(eid);
        eid = next_edge[get_ordered_vertex(eid, 1)];
        assert(eid != INVALID);
    }
    assert(eid == start_edge);

    // Update corner boundary vertices.
    for (size_t i = 0; i < num_bd_edges; i++) {
        const auto eid = combined_face.edges[i];
        const auto vid = get_ordered_vertex(eid, 0);

        auto& v = vertices[vid];
        if (v[0] <= 2 && v[1] <= 2) {
            v[2] = material_index;
        } else if (v[0] <= 2 && v[2] <= 2) {
            v[1] = material_index;
        } else if (v[1] <= 2 && v[2] <= 2) {
            v[0] = material_index;
        }
    }

    // Remove collinear boundary edges.
    auto vertex_on_boundary = [&](const Joint<2>& v) {
        size_t num_bd_material = 0;
        if (v[0] <= 2) num_bd_material++;
        if (v[1] <= 2) num_bd_material++;
        if (v[2] <= 2) num_bd_material++;
        return num_bd_material == 1;
    };

    bool skip_next = false;
    size_t output_edge_count = 0;
    for (size_t i = 0; i < num_bd_edges; i++) {
        if (skip_next) {
            skip_next = false;
            continue;
        }
        const auto eid = combined_face.edges[i];
        const auto v1 = get_ordered_vertex(eid, 1);

        const auto& v = vertices[v1];
        if (v[0] != material_index && v[1] != material_index && v[2] != material_index) {
            assert(vertex_on_boundary(v));
            logger().debug("remove boundary vertex {}: {} {} {}", v1, v[0], v[1], v[2]);

            const auto next_eid = combined_face.edges[(i + 1) % num_bd_edges];
            const auto v0 = get_ordered_vertex(eid, 0);
            const auto v2 = get_ordered_vertex(next_eid, 1);

            auto& e = edges[eid];
            e.vertices = {v0, v2};
            if (e.negative_material_label != material_index) {
                std::swap(e.positive_material_label, e.negative_material_label);
            }
            skip_next = true;
        }
        if (i != output_edge_count) {
            combined_face.edges[output_edge_count] = combined_face.edges[i];
        }
        output_edge_count++;
    }
    combined_face.edges.resize(output_edge_count);
    faces.push_back(std::move(combined_face));
    return faces.size() - 1;
}

} // namespace simplicial_arrangement
