#include "mi_cut_3_face.h"
#include "utils.h"

#include <absl/container/flat_hash_map.h>

#include <algorithm>
#include <array>
#include <vector>

namespace simplicial_arrangement {

std::array<size_t, 3> mi_cut_3_face(MIComplex<3>& mi_complex,
    size_t cid,
    size_t material_index,
    const std::vector<std::array<size_t, 3>>& subfaces)
{
    auto& edges = mi_complex.edges;
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

    auto compute_cut_edge_orientation = [&](size_t fid,
                                            const std::array<size_t, 3>& subface) -> bool {
        assert(subface[2] != INVALID);
        const auto& f = faces[fid];
        bool s = c.material_label == f.positive_material_label;

        if (subface[0] == INVALID || subface[1] == INVALID) {
            // Intersection edge is on the boundary of the face.
            auto itr = std::find(f.edges.begin(), f.edges.end(), subface[2]);
            assert(itr != f.edges.end());
            size_t curr_i = itr - f.edges.begin();
            size_t next_i = (curr_i + 1) % f.edges.size();

            const auto& curr_e = edges[f.edges[curr_i]];
            const auto& next_e = edges[f.edges[next_i]];
            bool edge_is_consistent_with_face = (curr_e.vertices[1] == next_e.vertices[0] ||
                                                 curr_e.vertices[1] == next_e.vertices[1]);

            bool on_positive_side = subface[0] != INVALID;
            logger().debug("cell/face: {}, face/edge: {}, face/cut plane: {}",
                s,
                edge_is_consistent_with_face,
                on_positive_side);

            uint8_t key = 0;
            if (s) key++;
            if (edge_is_consistent_with_face) key++;
            if (on_positive_side) key++;
            return key % 2 != 0;
        } else {
            // Intersection edge is a cross cut.
            return s;
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
            cut_edge_orientations.push_back(compute_cut_edge_orientation(fid, subface));
        }
    }

    if (positive_subfaces.empty() && negative_subfaces.empty()) {
        // The added material and the current material are identical over the
        // whole cell.
        return {INVALID, INVALID, INVALID};
    } else if (positive_subfaces.empty()) {
        cells[cid].material_label = material_index;
        return {INVALID, cid, cut_face_id};
    } else if (negative_subfaces.empty()) {
        return {cid, INVALID, cut_face_id};
    }

    // Chian cut edges into a loop
    {
        size_t num_cut_edges = cut_edges.size();
        assert(num_cut_edges >= 3);
        absl::flat_hash_map<size_t, size_t> v2e;
        v2e.reserve(num_cut_edges);
        for (size_t i = 0; i < num_cut_edges; i++) {
            const auto eid = cut_edges[i];
            const auto& e = edges[eid];
            if (cut_edge_orientations[i]) {
                v2e[e.vertices[0]] = i;
            } else {
                v2e[e.vertices[1]] = i;
            }
        }
        std::vector<size_t> chained_cut_edges;
        chained_cut_edges.reserve(num_cut_edges);
        chained_cut_edges.push_back(0);
        while (chained_cut_edges.size() < num_cut_edges) {
            const size_t i = chained_cut_edges.back();
            const auto& e = edges[cut_edges[i]];
            const size_t vid = cut_edge_orientations[i] ? e.vertices[1] : e.vertices[0];
            const auto itr = v2e.find(vid);
            assert(itr != v2e.end());
            const size_t next_i = itr->second;
            if (cut_edges[next_i] == cut_edges[chained_cut_edges.front()]) {
                break;
            }
            chained_cut_edges.push_back(next_i);
        }
        std::transform(chained_cut_edges.begin(),
            chained_cut_edges.end(),
            chained_cut_edges.begin(),
            [&](size_t i) { return cut_edges[i]; });
        std::swap(cut_edges, chained_cut_edges);
    }

    // Cross cut.
    assert(!cut_edges.empty());
    MIFace<3> cut_face;
    cut_face.edges = std::move(cut_edges);
    cut_face.positive_material_label = material_index;
    cut_face.negative_material_label = c.material_label;
    faces.push_back(std::move(cut_face));
    cut_face_id = faces.size() - 1;
    logger().debug("Adding cut face: {}", cut_face_id);
    assert(utils::edges_are_ordered(edges, faces[cut_face_id].edges));

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

    cells.push_back(std::move(positive_cell));
    cells.push_back(std::move(negative_cell));
    logger().debug("Adding positive subcell: {}", cells.size() - 2);
    logger().debug("Adding negative subcell: {}", cells.size() - 1);
    return {cells.size() - 2, cells.size() - 1, cut_face_id};
}

} // namespace simplicial_arrangement
