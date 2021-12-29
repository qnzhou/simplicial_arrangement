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
    const size_t num_vertices = mi_complex.vertices.size();
    const size_t num_edges = mi_complex.edges.size();
    const size_t num_faces = mi_complex.faces.size();
    const size_t num_cells = mi_complex.cells.size();
    logger().info("adding material {}", material_index);
    logger().debug("Before: {} {} {} {}", num_vertices, num_edges, num_faces, num_cells);

    auto& vertices = mi_complex.vertices;
    auto& edges = mi_complex.edges;
    auto& faces = mi_complex.faces;
    auto& cells = mi_complex.cells;

    // Reserve capacity. TODO: check if this helps.
    vertices.reserve(num_vertices + num_edges);
    edges.reserve(num_edges * 2);
    faces.reserve(num_faces * 2);
    cells.reserve(num_cells * 2);

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

    // Step 4: handle 3-faces.
    std::vector<std::array<size_t, 3>> subcells;
    subcells.reserve(num_cells);
    for (size_t i = 0; i < num_cells; i++) {
        subcells.push_back(mi_cut_3_face(builder, mi_complex, i, material_index, subfaces));
    }

    // Step 5: combine negative cells into a single cell.
    //{
    //    size_t combined_negative_cell =
    //        mi_union_negative_faces(builder, mi_complex, material_index, subfaces);
    //    std::vector<bool> to_keep(faces.size(), false);
    //    if (combined_negative_cell != INVALID) {
    //        to_keep[combined_negative_cell] = true;
    //    }
    //    std::for_each(subfaces.begin(), subfaces.end(), [&](const auto& subface) {
    //        if (subface[0] != INVALID) {
    //            to_keep[subface[0]] = true;
    //        }
    //    });
    //    shrink(faces, [&](size_t i) { return to_keep[i]; });
    //}

    // Step 6: consolidate.
    consolidate(mi_complex);
    logger().debug("After: {} {} {}",
        mi_complex.vertices.size(),
        mi_complex.edges.size(),
        mi_complex.faces.size());
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
