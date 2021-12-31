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
    auto& faces = mi_complex.faces;
    auto& cells = mi_complex.cells;

    const auto& c = cells[cid];
    const size_t num_bd_faces = c.faces.size();
    logger().debug("cell {} has {} bd faces", cid, num_bd_faces);

    size_t cut_face_id = INVALID;
    std::vector<size_t> positive_subfaces;
    std::vector<size_t> negative_subfaces;
    std::vector<size_t> cut_edges;
    positive_subfaces.reserve(num_bd_faces);
    negative_subfaces.reserve(num_bd_faces);
    cut_edges.reserve(num_bd_faces);

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
        }
    }

    if (positive_subfaces.empty()) {
        return {INVALID, cid, cut_face_id};
    } else if (negative_subfaces.empty()) {
        return {cid, INVALID, cut_face_id};
    }

    // Cross cut.
    assert(!cut_edges.empty());
    MIFace<3> cut_face;
    cut_face.edges = std::move(cut_edges);
    cut_face.positive_material_label = material_index;
    cut_face.negative_material_label = c.material_label;
    faces.push_back(std::move(cut_face));
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
    return {cells.size() - 2, cells.size() - 1, cut_face_id};
}

} // namespace simplicial_arrangement
