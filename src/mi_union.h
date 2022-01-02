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
template <int DIM, typename EdgeMap = absl::flat_hash_map<size_t, size_t>>
size_t mi_union_1_faces(MIComplex<DIM>& mi_complex,
    size_t material_index,
    const std::vector<size_t>& bd_edge_indices,
    EdgeMap* edge_map = nullptr)
{
    const size_t num_bd_edges = bd_edge_indices.size();
    if (num_bd_edges == 0) return INVALID;

    [[maybe_unused]] size_t edge_group_key = INVALID;
    if constexpr (DIM == 3) {
        if (edge_map != nullptr) {
            // Check if this set of edges has already been unioned.
            edge_group_key = *std::min_element(bd_edge_indices.begin(), bd_edge_indices.end());
            auto itr = edge_map->find(edge_group_key);
            if (itr != edge_map->end()) {
                logger().debug("Edge group with key {} has been unioned already as edge {}.",
                    edge_group_key,
                    itr->second);
                return itr->second;
            }
        }
    }

    auto& vertices = mi_complex.vertices;
    auto& edges = mi_complex.edges;

    auto compute_bd_vertices = [&]() {
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
        return bd_vertices;
    };

    size_t combined_eid = INVALID;
    if (num_bd_edges > 1) {
        MIEdge<DIM> combined_edge;
        auto bd_vertices = compute_bd_vertices();
        size_t counter = 0;
        for (auto vid : bd_vertices) {
            combined_edge.vertices[counter] = vid;
            counter++;
        }
        edges.push_back(std::move(combined_edge));
        combined_eid = edges.size() - 1;
        logger().debug("Adding combined edge: {}", combined_eid);
    } else {
        assert(bd_edge_indices.size() == 1);
        combined_eid = bd_edge_indices[0];
    }
    assert(combined_eid != INVALID);
    auto& combined_edge = edges[combined_eid];

    if constexpr (DIM == 3) {
        if (edge_map != nullptr) {
            // Register combined edge in edge map.
            assert(edge_group_key != INVALID);
            edge_map->insert({edge_group_key, combined_eid});
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

    return combined_eid;
}

template <int DIM, typename EdgeMap = absl::flat_hash_map<size_t, size_t>>
size_t mi_union_2_faces(MIComplex<DIM>& mi_complex,
    size_t material_index,
    const std::vector<size_t>& bd_face_indices,
    EdgeMap* edge_map = nullptr)
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
        size_t out_eid = mi_union_1_faces(mi_complex, material_index, edge_group, edge_map);
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

    absl::flat_hash_map<size_t, size_t> edge_map;
    edge_map.reserve(4);

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
        size_t combined_fid = mi_union_2_faces(mi_complex, material_index, face_group, &edge_map);
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

} // namespace simplicial_arrangement
