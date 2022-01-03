#include "ar_cut_3_face.h"

#include <absl/container/flat_hash_map.h>

#include <algorithm>
#include <array>
#include <vector>

namespace simplicial_arrangement {

std::array<size_t, 3> ar_cut_3_face(ARComplex<3>& ar_complex,
    size_t cid,
    size_t plane_index,
    const std::vector<std::array<size_t, 3>>& subfaces)
{
    auto& faces = ar_complex.faces;
    auto& cells = ar_complex.cells;

    const auto& c = cells[cid];
    const size_t num_bd_faces = c.faces.size();
    logger().debug("cell {} has {} bd faces", cid, num_bd_faces);

    size_t cut_face_id = INVALID;
    std::vector<size_t> positive_subfaces;
    std::vector<size_t> negative_subfaces;
    std::vector<size_t> cut_edges;
    positive_subfaces.reserve(num_bd_faces + 1);
    negative_subfaces.reserve(num_bd_faces + 1);
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

    if (positive_subfaces.empty() && negative_subfaces.empty()) {
        // The implicit function is identical over the whole cell.
        return {INVALID, INVALID, INVALID};
    } else if (positive_subfaces.empty()) {
        cells[cid].signs[plane_index] = false;
        return {INVALID, cid, cut_face_id};
    } else if (negative_subfaces.empty()) {
        cells[cid].signs[plane_index] = true;
        return {cid, INVALID, cut_face_id};
    }

    // Chain cut edges into a loop.
    {
        const auto& edges = ar_complex.edges;
        size_t num_cut_edges = cut_edges.size();
        assert(num_cut_edges >= 3);
        size_t num_vertices = ar_complex.vertices.size();
        std::vector<size_t> v2e(num_vertices, INVALID);
        for (size_t i=0; i<num_cut_edges; i++) {
            const auto eid = cut_edges[i];
            const auto& e = edges[eid];
            v2e[e.vertices[0]] = eid;
        }
        std::vector<size_t> chained_cut_edges;
        chained_cut_edges.reserve(num_cut_edges);
        chained_cut_edges.push_back(cut_edges.front());
        while(chained_cut_edges.size() < num_cut_edges) {
            const auto& e = edges[chained_cut_edges.back()];
            const size_t next_e = v2e[e.vertices[1]];
            assert(next_e != INVALID);
            if (next_e == chained_cut_edges.front()) {
                break;
            }
            chained_cut_edges.push_back(next_e);
        }
        std::swap(cut_edges, chained_cut_edges);
    }

    // Cross cut.
    assert(!cut_edges.empty());
    ARFace<3> cut_face;
    cut_face.edges = std::move(cut_edges);
    cut_face.supporting_plane = plane_index;
    faces.push_back(std::move(cut_face));
    cut_face_id = faces.size() - 1;
    logger().debug("Adding cut face: {}", cut_face_id);

    // Generate positive and negative subcell.
    ARCell<3> positive_cell, negative_cell;
    positive_cell.faces.reserve(positive_subfaces.size() + 1);
    negative_cell.faces.reserve(negative_subfaces.size() + 1);

    positive_subfaces.push_back(cut_face_id);
    positive_cell.faces = std::move(positive_subfaces);
    positive_cell.signs = c.signs;
    positive_cell.signs[plane_index] = true;

    negative_subfaces.push_back(cut_face_id);
    negative_cell.faces = std::move(negative_subfaces);
    negative_cell.signs = c.signs;
    negative_cell.signs[plane_index] = false;

    cells.push_back(std::move(positive_cell));
    cells.push_back(std::move(negative_cell));
    logger().debug("Adding positive subcell: {}", cells.size() - 2);
    logger().debug("Adding negative subcell: {}", cells.size() - 1);

    {
        assert(cut_face_id != INVALID);
        auto& cut_f = faces[cut_face_id];
        cut_f.positive_cell = cells.size() - 2;
        cut_f.negative_cell = cells.size() - 1;
    }
    return {cells.size() - 2, cells.size() - 1, cut_face_id};
}

} // namespace simplicial_arrangement
