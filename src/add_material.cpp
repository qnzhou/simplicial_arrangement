#include "add_material.h"
#include "MaterialInterfaceBuilder.h"
#include "mi_cut.h"
#include "mi_union.h"
#include "utils.h"

#include <implicit_predicates/implicit_predicates.h>

#include <optional>

namespace {
using namespace simplicial_arrangement;

/**
 * Remove unused verices and faces.
 */
template <int DIM>
void consolidate(MIComplex<DIM>& mi_complex)
{
    // Shrink faces.
    if constexpr (DIM == 3) {
        std::vector<bool> active_faces(mi_complex.faces.size(), false);
        for (auto& c : mi_complex.cells) {
            for (auto fid : c.faces) {
                active_faces[fid] = true;
            }
        }

        auto index_map = shrink(mi_complex.faces, [&](size_t fid) { return active_faces[fid]; });

        for (auto& c : mi_complex.cells) {
            std::transform(c.faces.begin(), c.faces.end(), c.faces.begin(), [&](size_t i) {
                assert(index_map[i] != INVALID);
                return index_map[i];
            });
        }
    }

    // Shrink edges.
    {
        std::vector<bool> active_edges(mi_complex.edges.size(), false);
        for (auto& f : mi_complex.faces) {
            for (auto eid : f.edges) {
                active_edges[eid] = true;
            }
        }

        auto index_map = shrink(mi_complex.edges, [&](size_t eid) { return active_edges[eid]; });

        for (auto& f : mi_complex.faces) {
            std::transform(f.edges.begin(), f.edges.end(), f.edges.begin(), [&](size_t i) {
                assert(index_map[i] != INVALID);
                return index_map[i];
            });
        }
    }

    // Shrink vertices.
    {
        std::vector<bool> active_vertices(mi_complex.vertices.size(), false);
        for (auto& e : mi_complex.edges) {
            for (auto vid : e.vertices) {
                assert(vid != INVALID);
                active_vertices[vid] = true;
            }
        }

        auto index_map =
            shrink(mi_complex.vertices, [&](size_t vid) { return active_vertices[vid]; });

        for (auto& e : mi_complex.edges) {
            e.vertices[0] = index_map[e.vertices[0]];
            e.vertices[1] = index_map[e.vertices[1]];
        }
    }
}

template <typename Scalar>
void add_material(
    MaterialInterfaceBuilder<Scalar, 2>& builder, MIComplex<2>& mi_complex, size_t material_index)
{
    const size_t num_vertices = mi_complex.vertices.size();
    const size_t num_edges = mi_complex.edges.size();
    const size_t num_faces = mi_complex.faces.size();
    logger().info("adding material {}", material_index);
    logger().debug("Before: {} {} {}", num_vertices, num_edges, num_faces);

    auto& vertices = mi_complex.vertices;
    auto& edges = mi_complex.edges;
    auto& faces = mi_complex.faces;

    // Reserve capacity.
    vertices.reserve(num_vertices + num_edges);
    edges.reserve(num_edges * 2);
    faces.reserve(num_faces * 2);

    // Step 1: handle 0-faces.
    std::vector<int8_t> orientations;
    orientations.reserve(num_vertices);
    for (size_t i = 0; i < num_vertices; i++) {
        orientations.push_back(mi_cut_0_face(builder, mi_complex, i, material_index));
    }

    // Step 2: handle 1-faces.
    std::vector<std::array<size_t, 3>> subedges;
    subedges.reserve(num_edges);
    for (size_t i = 0; i < num_edges; i++) {
        subedges.push_back(mi_cut_1_face(builder, mi_complex, i, material_index, orientations));
    }

    // Step 3: handle 2-faces.
    std::vector<std::array<size_t, 3>> subfaces;
    subfaces.reserve(num_faces);
    for (size_t i = 0; i < num_faces; i++) {
        subfaces.push_back(
            mi_cut_2_face(builder, mi_complex, i, material_index, orientations, subedges));
    }

    // Step 4: combine negative cells into a single cell.
    {
        size_t combined_negative_cell =
            mi_union_negative_faces(builder, mi_complex, material_index, subfaces);
        std::vector<bool> to_keep(faces.size(), false);
        if (combined_negative_cell != INVALID) {
            to_keep[combined_negative_cell] = true;
        }
        std::for_each(subfaces.begin(), subfaces.end(), [&](const auto& subface) {
            if (subface[0] != INVALID) {
                to_keep[subface[0]] = true;
            }
        });
        shrink(faces, [&](size_t i) { return to_keep[i]; });
    }

    // Step 5: consolidate.
    consolidate(mi_complex);
    logger().debug("After: {} {} {}",
        mi_complex.vertices.size(),
        mi_complex.edges.size(),
        mi_complex.faces.size());
}

template <typename Scalar>
void add_material(
    MaterialInterfaceBuilder<Scalar, 3>& builder, MIComplex<3>& mi_complex, size_t material_index)
{
    // auto& material_interface = builder.get_material_interface();
    // const auto& material = builder.get_material(material_index);
    // const size_t num_vertices = material_interface.vertices.size();
    // const size_t num_faces = material_interface.faces.size();

    // auto& vertices = material_interface.vertices;
    // auto& edges = material_interface.edges;
    // auto& faces = material_interface.faces;

    //// Step 1: cut 0-face.
    // std::vector<implicit_predicates::Orientation> orientations;
    // orientations.reserve(num_vertices);
    // for (size_t i = 0; i < num_vertices; i++) {
    //    const auto& p = vertices[i];
    //    orientations.push_back(compute_orientation(builder, p, material));
    //    assert(orientations.back() != implicit_predicates::INVALID);
    //    logger().debug("vertex {}: {}, {}, {}", i, p[0], p[1], p[2]);
    //    logger().debug("orientation: {}", orientations.back());
    //}

    //// Step 2: cut 2-face.
    // std::vector<size_t> positive_subfaces;
    // std::vector<size_t> negative_subfaces;
    // std::vector<std::array<size_t, 2>> cut_edges;
    // positive_subfaces.reserve(num_faces);
    // negative_subfaces.reserve(num_faces);
    // cut_edges.reserve(num_faces);

    // auto cut_face = [&](size_t i) {
    //    const auto& f = faces[i];

    //    // Check for easy cases first.
    //    bool all_positive = true;
    //    bool all_negative = true;
    //    bool all_zero = true;
    //    for (auto eid : f.edges) {
    //        const auto& vids = edges[eid].vertices;
    //        if (orientations[vids[0]] > 0 && orientations[vids[1]] > 0) {
    //            all_negative = false;
    //            all_zero = false;
    //        } else if (orientations[vids[0]] < 0 && orientations[vids[1]] < 0) {
    //            all_positive = false;
    //            all_zero = false;
    //        } else if (orientations[vids[0]] == 0 && orientations[vids[1]] == 0){
    //            all_positive = false;
    //            all_negative = false;
    //        } else {
    //            all_positive = false;
    //            all_negative = false;
    //            all_zero = false;
    //            break;
    //        }
    //    }
    //    if (all_positive) {
    //        positive_subfaces.push_back(i);
    //        negative_subfaces.push_back(INVALID);
    //        cut_edges.push_back({INVALID, INVALID});
    //        return;
    //    } else if (all_negative) {
    //        positive_subfaces.push_back(INVALID);
    //        negative_subfaces.push_back(i);
    //        cut_edges.push_back({INVALID, INVALID});
    //        return;
    //    } else if (all_zero) {
    //        positive_subfaces.push_back(INVALID);
    //        negative_subfaces.push_back(INVALID);
    //        cut_edges.push_back({INVALID, INVALID});
    //        return;
    //    }

    //    // The cut is either a cross cut or is tangent at 1 or 2 vertices.
    //    //const size_t num_bd_edges = f.edges.size();
    //    //size_t first_positive = INVALID;
    //    ////size_t first_negative = INVALID;
    //    //std::array<size_t, 2> intersection{INVALID, INVALID};
    //    //for (size_t j = 0; j < num_bd_vertices; j++) {
    //    //    const size_t v0 = f.vertices[j];
    //    //    const size_t v1 = f.vertices[(j + 1) % num_bd_vertices];
    //    //    const auto o0 = orientations[v0];
    //    //    const auto o1 = orientations[v1];
    //    //    if (o0 == 0 && o1 == 0) {
    //    //        intersection = {v0, v1};
    //    //        continue;
    //    //    } else if (o0 <= 0 && o1 > 0) {
    //    //        first_positive = j;
    //    //        // TODO.
    //    //    }
    //    //}
    //};

    // for (size_t i = 0; i < num_faces; i++) {
    //    cut_face(i);
    //}
}

} // namespace

namespace simplicial_arrangement::internal {

void add_material(
    MaterialInterfaceBuilder<double, 2>& builder, MIComplex<2>& mi_complex, size_t material_index)
{
    ::add_material(builder, mi_complex, material_index);
}

void add_material(
    MaterialInterfaceBuilder<Int, 2>& builder, MIComplex<2>& mi_complex, size_t material_index)
{
    ::add_material(builder, mi_complex, material_index);
}

void add_material(
    MaterialInterfaceBuilder<double, 3>& builder, MIComplex<3>& mi_complex, size_t material_index)
{
    ::add_material(builder, mi_complex, material_index);
}

void add_material(
    MaterialInterfaceBuilder<Int, 3>& builder, MIComplex<3>& mi_complex, size_t material_index)
{
    ::add_material(builder, mi_complex, material_index);
}

} // namespace simplicial_arrangement::internal
