#include "add_material.h"
#include "MaterialInterfaceBuilder.h"
#include "consolidate_complex.h"
#include "mi_cut.h"
#include "mi_union.h"

#include <array>
#include <vector>

namespace {
using namespace simplicial_arrangement;

template <typename Scalar>
void add_material(
    MaterialInterfaceBuilder<Scalar, 2>& builder, MIComplex<2>& mi_complex, size_t material_index)
{
    const size_t num_vertices = mi_complex.vertices.size();
    const size_t num_edges = mi_complex.edges.size();
    const size_t num_faces = mi_complex.faces.size();
    logger().debug("adding material {}", material_index);
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
        orientations.push_back(
            mi_cut_0_face(builder.get_material_repo(), mi_complex, i, material_index));
    }

    // Step 2: handle 1-faces.
    std::vector<std::array<size_t, 3>> subedges;
    subedges.reserve(num_edges);
    for (size_t i = 0; i < num_edges; i++) {
        subedges.push_back(mi_cut_1_face(mi_complex, i, material_index, orientations));
    }

    // Step 3: handle 2-faces.
    std::vector<std::array<size_t, 3>> subfaces;
    subfaces.reserve(num_faces);
    for (size_t i = 0; i < num_faces; i++) {
        subfaces.push_back(mi_cut_2_face(mi_complex, i, material_index, orientations, subedges));
    }

    // Step 4: combine negative cells into a single cell.
    {
        std::vector<size_t> negative_subfaces;
        negative_subfaces.reserve(subfaces.size());
        for (const auto& subface : subfaces) {
            if (subface[1] != INVALID) {
                negative_subfaces.push_back(subface[1]);
            }
        }
        size_t combined_negative_cell =
            mi_union_2_faces(mi_complex, material_index, negative_subfaces);
        std::vector<bool> to_keep(faces.size(), false);
        if (combined_negative_cell != INVALID) {
            to_keep[combined_negative_cell] = true;
        }
        std::for_each(subfaces.begin(), subfaces.end(), [&](const auto& subface) {
            if (subface[0] != INVALID) {
                to_keep[subface[0]] = true;
            }
        });
        utils::shrink(faces, [&](size_t i) { return to_keep[i]; });
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
    logger().debug("adding material {}", material_index);
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
        orientations.push_back(
            mi_cut_0_face(builder.get_material_repo(), mi_complex, i, material_index));
    }

    // Step 2: handle 1-faces.
    std::vector<std::array<size_t, 3>> subedges;
    subedges.reserve(num_edges);
    for (size_t i = 0; i < num_edges; i++) {
        subedges.push_back(mi_cut_1_face(mi_complex, i, material_index, orientations));
    }

    // Step 3: handle 2-faces.
    std::vector<std::array<size_t, 3>> subfaces;
    subfaces.reserve(num_faces);
    for (size_t i = 0; i < num_faces; i++) {
        subfaces.push_back(mi_cut_2_face(mi_complex, i, material_index, orientations, subedges));
    }

    // Step 4: handle 3-faces.
    std::vector<std::array<size_t, 3>> subcells;
    subcells.reserve(num_cells);
    for (size_t i = 0; i < num_cells; i++) {
        subcells.push_back(mi_cut_3_face(mi_complex, i, material_index, subfaces));
    }

    // Step 5: combine negative cells into a single cell.
    {
        std::vector<size_t> negative_cells;
        negative_cells.reserve(num_cells);
        for (const auto& subcell : subcells) {
            if (subcell[1] == INVALID) continue;
            negative_cells.push_back(subcell[1]);
        }
        size_t combined_negative_cell =
            mi_union_3_faces(mi_complex, material_index, negative_cells);
        std::vector<bool> to_keep(cells.size(), false);
        if (combined_negative_cell != INVALID) {
            to_keep[combined_negative_cell] = true;
        }
        std::for_each(subcells.begin(), subcells.end(), [&](const auto& subcell) {
            if (subcell[0] != INVALID) {
                to_keep[subcell[0]] = true;
            }
        });
        utils::shrink(cells, [&](size_t i) { return to_keep[i]; });
    }

    // Step 6: consolidate.
    consolidate(mi_complex);
    logger().debug("After: {} {} {} {}",
        mi_complex.vertices.size(),
        mi_complex.edges.size(),
        mi_complex.faces.size(),
        mi_complex.cells.size());
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
