#pragma once

#include "MIComplex.h"
#include "MaterialInterfaceBuilder.h"

#include <absl/container/flat_hash_map.h>
#include <absl/container/flat_hash_set.h>
#include <array>
#include <vector>

namespace simplicial_arrangement {

/**
 * Combined multiple edge into a single edge.  This operation only makes sense
 * for boundary edges.  This function also update the boundary vertex and edge's
 * supporting materials.
 */
template <int DIM>
size_t mi_union_1_faces(
    MIComplex<DIM>& mi_complex, size_t material_index, const std::vector<size_t>& bd_edge_indices)
{
    const size_t num_bd_edges = bd_edge_indices.size();
    if (num_bd_edges == 0) return INVALID;

    auto& vertices = mi_complex.vertices;
    auto& edges = mi_complex.edges;
    absl::flat_hash_set<size_t> bd_vertices;
    bd_vertices.reserve(num_bd_edges + 1);

    auto toggle = [&](size_t vid) {
        auto itr = bd_vertices.find(vid);
        if (itr == bd_vertices.end()) {
            bd_vertices.insert(vid);
        } else {
            bd_vertices.erase(itr);
        }
    };

    for (size_t eid : bd_edge_indices) {
        const auto& e = edges[eid];
        toggle(e.vertices[0]);
        toggle(e.vertices[1]);
    }
    assert(bd_vertices.size() == 2);

    MIEdge<DIM> combined_edge;
    {
        size_t counter = 0;
        for (auto vid : bd_vertices) {
            combined_edge.vertices[counter] = vid;
            counter++;
        }
    }

    // Correct vertex material label.
    auto correct_vertex_material_label = [&](size_t vid) {
        auto& v = vertices[vid];
        if constexpr (DIM == 2) {
            if (v[0] <= DIM && v[1] <= DIM) {
                v[2] = material_index;
            } else if (v[0] <= DIM && v[2] <= DIM) {
                v[1] = material_index;
            } else if (v[1] <= DIM && v[2] <= DIM) {
                v[0] = material_index;
            }
        } else {
            if (v[0] <= DIM && v[1] <= DIM && v[2] <= DIM) {
                v[3] = material_index;
            } else if (v[0] <= DIM && v[1] <= DIM && v[3] <= DIM) {
                v[2] = material_index;
            } else if (v[0] <= DIM && v[2] <= DIM && v[3] <= DIM) {
                v[1] = material_index;
            } else if (v[1] <= DIM && v[2] <= DIM && v[3] <= DIM) {
                v[0] = material_index;
            }
        }
    };
    correct_vertex_material_label(combined_edge.vertices[0]);
    correct_vertex_material_label(combined_edge.vertices[1]);

    if constexpr (DIM == 3) {
        combined_edge.supporting_materials = edges[bd_edge_indices.front()].supporting_materials;
        auto correct_edge_material = [&](size_t i) {
            if (combined_edge.supporting_materials[i] > 3) {
                combined_edge.supporting_materials[i] = material_index;
                assert(combined_edge.supporting_materials[(i + 1) % 3] <= 3);
                assert(combined_edge.supporting_materials[(i + 2) % 3] <= 3);
            }
        };
        correct_edge_material(0);
        correct_edge_material(1);
        correct_edge_material(2);
    } else {
        combined_edge.positive_material_label =
            edges[bd_edge_indices.front()].positive_material_label;
        combined_edge.negative_material_label =
            edges[bd_edge_indices.front()].negative_material_label;
        if (combined_edge.positive_material_label > 2) {
            combined_edge.positive_material_label = material_index;
            assert(combined_edge.negative_material_label <= 2);
        } else if (combined_edge.negative_material_label > 2) {
            combined_edge.negative_material_label = material_index;
            assert(combined_edge.positive_material_label <= 2);
        } else {
            logger().error(
                "Both positive and negative materials of an edge are boundary material!");
            assert(false);
        }
    }

    edges.push_back(std::move(combined_edge));
    return edges.size() - 1;
}

template <int DIM>
size_t mi_union_2_faces(
    MIComplex<DIM>& mi_complex, size_t material_index, const std::vector<size_t>& bd_face_indices)
{
    const size_t num_bd_faces = bd_face_indices.size();
    if (num_bd_faces == 0) return INVALID;

    auto& edges = mi_complex.edges;
    auto& faces = mi_complex.faces;

    absl::flat_hash_set<size_t> bd_edges;
    size_t num_active_half_edges = 0;
    for (auto fid : bd_face_indices) {
        num_active_half_edges += faces[fid].edges.size();
    }
    bd_edges.reserve(num_active_half_edges);

    auto toggle = [&](size_t eid) {
        auto itr = bd_edges.find(eid);
        if (itr == bd_edges.end()) {
            bd_edges.insert(eid);
        } else {
            bd_edges.erase(itr);
        }
    };

    for (size_t fid : bd_face_indices) {
        const auto& f = faces[fid];
        for (size_t eid : f.edges) {
            toggle(eid);
        }

        if constexpr (DIM == 2) {
            assert(f.material_label == material_index);
        }
    }

    std::array<std::vector<size_t>, DIM*(DIM + 1) / 2> to_merge;
    std::for_each(to_merge.begin(), to_merge.end(), [&](std::vector<size_t>& entry) {
        entry.reserve(bd_edges.size());
    });

    auto sort3 = [](std::array<size_t, 3>& a) {
        if (a[0] > a[1]) std::swap(a[0], a[1]);
        if (a[0] > a[2]) std::swap(a[0], a[2]);
        if (a[1] > a[2]) std::swap(a[1], a[2]);
    };

    for (size_t eid : bd_edges) {
        auto& e = edges[eid];
        if constexpr (DIM == 2) {
            if (e.positive_material_label <= 2) {
                assert(e.negative_material_label > 2);
                to_merge[e.positive_material_label].push_back(eid);
            } else if (e.negative_material_label <= 2) {
                assert(e.positive_material_label > 2);
                to_merge[e.negative_material_label].push_back(eid);
            }
        } else {
            sort3(e.supporting_materials);
            size_t m0 = e.supporting_materials[0];
            size_t m1 = e.supporting_materials[1];
            if (m1 <= 3) {
                assert(e.supporting_materials[2] > 3);
                if (m0 == 0 && m1 == 1) {
                    to_merge[0].push_back(eid);
                } else if (m0 == 0 && m1 == 2) {
                    to_merge[1].push_back(eid);
                } else if (m0 == 0 && m1 == 3) {
                    to_merge[2].push_back(eid);
                } else if (m0 == 1 && m1 == 2) {
                    to_merge[3].push_back(eid);
                } else if (m0 == 1 && m1 == 3) {
                    to_merge[4].push_back(eid);
                } else if (m0 == 2 && m1 == 3) {
                    to_merge[5].push_back(eid);
                }
            }
        }
    }
    for (const auto& edge_group : to_merge) {
        if (edge_group.empty()) continue;
        size_t out_eid = mi_union_1_faces(mi_complex, material_index, edge_group);
        assert(out_eid != INVALID);
        for (auto eid : edge_group) {
            bd_edges.erase(eid);
        }
        bd_edges.insert(out_eid);
    }

    size_t num_bd_edges = bd_edges.size();
    MIFace<DIM> combined_face;
    combined_face.edges.reserve(num_bd_edges);
    for (size_t eid : bd_edges) {
        combined_face.edges.push_back(eid);
    }
    if constexpr (DIM == 2) {
        combined_face.material_label = material_index;
    } else {
        // For 3D combined face must be on a tet boundary triangle.
        const auto& old_face = faces[bd_face_indices.front()];
        if (old_face.positive_material_label <= 3) {
            assert(old_face.negative_material_label > 3);
            combined_face.positive_material_label = old_face.positive_material_label;
            combined_face.negative_material_label = material_index;
        } else {
            assert(old_face.positive_material_label > 3);
            assert(old_face.negative_material_label <= 3);
            combined_face.positive_material_label = material_index;
            combined_face.negative_material_label = old_face.negative_material_label;
        }
    }

    faces.push_back(std::move(combined_face));
    return faces.size() - 1;
}

template <int DIM>
size_t mi_union_3_faces(
    MIComplex<DIM>& mi_complex, size_t material_index, const std::vector<size_t>& bd_cell_indices)
{
    static_assert(DIM == 3);
    const size_t num_bd_cells = bd_cell_indices.size();
    if (num_bd_cells == 0) return INVALID;

    auto& faces = mi_complex.faces;
    auto& cells = mi_complex.cells;

    // Extract a set of bounding faces to the combined cell.
    absl::flat_hash_set<size_t> bd_faces;
    size_t num_active_half_faces = 0;
    for (auto cid : bd_cell_indices) {
        num_active_half_faces += cells[cid].faces.size();
    }
    bd_faces.reserve(num_active_half_faces);

    auto toggle = [&](size_t fid) {
        auto itr = bd_faces.find(fid);
        if (itr == bd_faces.end()) {
            bd_faces.insert(fid);
        } else {
            bd_faces.erase(itr);
        }
    };

    for (size_t cid : bd_cell_indices) {
        const auto& c = cells[cid];
        assert(c.material_label == material_index);
        for (size_t fid : c.faces) {
            toggle(fid);
        }
    }

    // Merge coplanar faces on the tet boundary.
    std::array<std::vector<size_t>, 4> to_merge;
    std::for_each(to_merge.begin(), to_merge.end(), [&](std::vector<size_t>& entry) {
        entry.reserve(bd_faces.size());
    });

    for (auto fid : bd_faces) {
        const auto& f = faces[fid];
        if (f.positive_material_label <= DIM) {
            to_merge[f.positive_material_label].push_back(fid);
        } else if (f.negative_material_label <= DIM) {
            to_merge[f.negative_material_label].push_back(fid);
        }
    }

    for (const auto& face_group : to_merge) {
        if (face_group.empty()) continue;
        size_t combined_fid = mi_union_2_faces(mi_complex, material_index, face_group);
        for (auto fid : face_group) {
            bd_faces.erase(fid);
        }
        bd_faces.insert(combined_fid);
    }

    // Generate combined cell.
    const size_t num_bd_faces = bd_faces.size();
    MICell<DIM> combined_cell;
    combined_cell.material_label = material_index;
    combined_cell.faces.reserve(num_bd_faces);
    for (auto fid : bd_faces) {
        combined_cell.faces.push_back(fid);
    }

    cells.push_back(std::move(combined_cell));
    return cells.size() - 1;
}

template <typename Scalar>
size_t mi_union_negative_faces(MaterialInterfaceBuilder<Scalar, 2>& builder,
    MIComplex<2>& mi_complex,
    size_t material_index,
    const std::vector<std::array<size_t, 3>>& subfaces)
{
    const size_t num_negative_faces = std::count_if(
        subfaces.begin(), subfaces.end(), [](const auto& val) { return val[1] != INVALID; });
    if (num_negative_faces == 0) return INVALID;
    logger().debug("combining {} negative 2D cells.", num_negative_faces);

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

    // Update corner boundary edges and vertices.
    for (size_t i = 0; i < num_bd_edges; i++) {
        const auto eid = combined_face.edges[i];
        const auto vid = get_ordered_vertex(eid, 0);

#ifndef NDEBUG
        auto& e = edges[eid];
        assert(e.positive_material_label == material_index ||
               e.negative_material_label == material_index);
#endif

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

// template <typename Scalar>
// size_t mi_union_negative_cells(MaterialInterfaceBuilder<Scalar, 3>& builder,
//    MIComplex<3>& mi_complex,
//    size_t material_index,
//    const std::vector<std::array<size_t, 3>>& subcells)
//{
//    const size_t num_negative_cells = std::count_if(
//        subcells.begin(), subcells.end(), [](const auto& val) { return val[1] != INVALID; });
//    if (num_negative_cells == 0) return INVALID;
//    logger().debug("combining {} negative 3D cells.", num_negative_cells);
//
//    auto& vertices = mi_complex.vertices;
//    auto& edges = mi_complex.edges;
//    auto& faces = mi_complex.faces;
//    auto& cells = mi_complex.cells;
//
//    MICell<3> combined_cell;
//    combined_cell.material_label = material_index;
//
//    const size_t total_num_vertices = mi_complex.vertices.size();
//    const size_t total_num_edges = mi_complex.edges.size();
//    const size_t total_num_faces = mi_complex.faces.size();
//
//    std::vector<bool> on_border(total_num_faces, false);
//
//    for (const auto& subcell : subcells) {
//        const auto cid = subcell[1]; // Negative subcell.
//        if (cid == INVALID) continue;
//        const auto& c = cells[cid];
//        for (const auto& fid : c.faces) {
//            on_border[fid] = !on_border[fid];
//        }
//    }
//
//    size_t num_bd_faces = std::count(on_border.begin(), on_border.end(), [](bool v) { return v;
//    }); std::vector<size_t> bd_faces; bd_faces.reserve(num_bd_faces); for (size_t i = 0; i <
//    total_num_faces; i++) {
//        if (!on_border[i]) continue;
//        bd_faces.push_back(i);
//    }
//
//    std::array<std::vector<size_t>, 4> bd_face_groups;
//    bd_face_groups[0].reserve(num_bd_faces);
//    bd_face_groups[1].reserve(num_bd_faces);
//    bd_face_groups[2].reserve(num_bd_faces);
//    bd_face_groups[3].reserve(num_bd_faces);
//    for (size_t fid : bd_faces) {
//        auto& f = faces[fid];
//        if (f.positive_material_label <= 3) {
//            f.negative_material_label = material_index;
//            bd_face_groups[f.positive_material_label].push_back(fid);
//        } else if (f.negative_material_label <= 3) {
//            f.positive_material_label = material_index;
//            bd_face_groups[f.negative_material_label].push_back(fid);
//        } else {
//            combined_cell.faces.push_back(fid);
//        }
//    }
//
//    for (size_t i = 0; i < 4; i++) {
//        if (bd_face_groups[i].empty()) continue;
//
//        // Correct for corner_vertices;
//        for (size_t fid : bd_face_groups[i]) {
//            // TODO.
//        }
//    }
//
//
//    return INVALID;
//}

} // namespace simplicial_arrangement
