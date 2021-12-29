#pragma once

#include "MIComplex.h"
#include "MaterialInterfaceBuilder.h"

#include <absl/container/flat_hash_map.h>
#include <array>
#include <vector>

namespace simplicial_arrangement {

template <typename Scalar>
std::array<size_t, 3> mi_cut_3_face(MaterialInterfaceBuilder<Scalar, 3>& builder,
    MIComplex<3>& mi_complex,
    size_t cid,
    size_t material_index,
    const std::vector<std::array<size_t, 3>>& subfaces)
{
    const auto& edges = mi_complex.edges;
    auto& faces = mi_complex.faces;
    auto& cells = mi_complex.cells;

    const auto& c = cells[cid];
    const size_t num_bd_faces = c.faces.size();
    logger().debug("cell {} has {} bd faces", cid, num_bd_faces);

    size_t cut_face_id = INVALID;
    std::vector<size_t> positive_subfaces;
    std::vector<size_t> negative_subfaces;
    std::vector<size_t> cut_edges;
    std::vector<bool> cut_edge_orientations;
    positive_subfaces.reserve(num_bd_faces);
    negative_subfaces.reserve(num_bd_faces);
    cut_edges.reserve(num_bd_faces);
    cut_edge_orientations.reserve(num_bd_faces);

    // Return true if edge is consistently orientated with respect to face.
    auto get_edge_orientation_with_respect_to_face = [&](size_t eid, size_t fid) {
        const auto& f = faces[fid];
        auto itr = std::find(f.edges.begin(), f.edges.end(), eid);
        assert(itr != f.edges.end());
        const size_t i = itr - f.edges.begin();
        const size_t j = (i + 1) % f.edges.size();
        const auto& ei = edges[eid];
        const auto& ej = edges[f.edges[j]];
        if (ei.vertices[0] == ej.vertices[0] || ei.vertices[0] == ej.vertices[1]) {
            return false;
        } else {
            assert(ei.vertices[1] == ej.vertices[0] || ei.vertices[1] == ej.vertices[1]);
            return true;
        }
    };

    for (auto fid : c.faces) {
        const auto& subface = subfaces[fid];
        if (subface[0] == INVALID && subface[1] == INVALID) {
            cut_face_id = fid;
        }
        if (subface[0] != INVALID) {
            positive_subfaces.push_back(subface[0]);
        }
        if (subface[1] != INVALID) {
            negative_subfaces.push_back(subface[1]);
        }
        if (subface[2] != INVALID) {
            cut_edges.push_back(subface[2]);
            if (subface[0] == INVALID) {
                // Face is on the negative side and cut edge is part of face boundary.
                cut_edge_orientations.push_back(
                    get_edge_orientation_with_respect_to_face(subface[0], fid));
            } else if (subface[1] == INVALID)
                // Face is on the positive side and cut edge is part of face boundary.
                cut_edge_orientations.push_back(
                    !get_edge_orientation_with_respect_to_face(subface[1], fid));
        } else {
            cut_edge_orientations.push_back(true);
        }
    }

    if (positive_subfaces.empty()) {
        return {INVALID, cid, cut_face_id};
    } else if (negative_subfaces.empty()) {
        return {cid, INVALID, cut_face_id};
    }

    // Cross cut.
    assert(!cut_edges.empty());
    const size_t num_cut_edges = cut_edges.size();
    MIFace<3> cut_face;
    {
        cut_face.edges.reserve(num_cut_edges);
        auto get_ordered_vertex = [&](size_t i, size_t j) {
            assert(i < num_cut_edges);
            assert(j < 2);
            const auto& e = edges[cut_edges[i]];
            if (cut_edge_orientations[i]) {
                return e.vertices[j];
            } else {
                return e.vertices[(1 + j) % 2];
            }
        };

        absl::flat_hash_map<size_t, size_t> next_edge;
        next_edge.reserve(num_cut_edges);
        for (size_t i = 0; i < num_cut_edges; i++) {
            next_edge[get_ordered_vertex(i, 0)] = i;
        }

        size_t id = 0;
        do {
            cut_face.edges.push_back(cut_edges[id]);
            id = next_edge[get_ordered_vertex(id, 1)];
            assert(cut_face.edges.size() < num_cut_edges);
        } while (id != 0);
        cut_face.positive_material_label = material_index;
        cut_face.negative_material_label = c.material_label;
    }
    faces.push_back(cut_face);
    cut_face_id = faces.size() - 1;

    // Generate positive and negative subcell.
    MICell<3> positive_cell, negative_cell;
    positive_cell.faces.reserve(positive_subfaces.size() + 1);
    negative_cell.faces.reserve(negative_subfaces.size() + 1);

    positive_cell.faces.insert(
        positive_cell.faces.begin(), positive_subfaces.begin(), positive_subfaces.end());
    positive_cell.faces.push_back(cut_face_id);
    negative_cell.faces.insert(
        negative_cell.faces.begin(), negative_subfaces.begin(), negative_subfaces.end());
    negative_cell.faces.push_back(cut_face_id);
    positive_cell.material_label = c.material_label;
    negative_cell.material_label = material_index;

    cells.push_back(positive_cell);
    cells.push_back(negative_cell);
    return {cells.size() - 1, cells.size() - 2, cut_face_id};
}

} // namespace simplicial_arrangement
