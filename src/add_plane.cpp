#include "add_plane.h"
#include "PlaneRepo.h"
#include "ar_cut.h"
#include "consolidate_complex.h"

#include <array>
#include <vector>

namespace {

template <typename Scalar>
void add_plane(const PlaneRepo<Scalar, 2>& repo, ARComplex<2>& ar_complex, size_t plane_index)
{
    const size_t num_vertices = ar_complex.vertices.size();
    const size_t num_edges = ar_complex.edges.size();
    const size_t num_faces = ar_complex.faces.size();
    logger().debug("adding material {}", plane_index);
    logger().debug("Before: {} {} {}", num_vertices, num_edges, num_faces);

    auto& vertices = ar_comlpex.vertices;
    auto& edges = ar_comlpex.edges;
    auto& faces = ar_comlpex.faces;

    // Reserve capacity.
    vertices.reserve(num_vertices + num_edges);
    edges.reserve(num_edges * 2);
    faces.reserve(num_faces * 2);

    // Step 1: handle 0-faces.
    std::vector<int8_t> orientations;
    orientations.reserve(num_vertices);
    for (size_t i = 0; i < num_vertices; i++) {
        orientations.push_back(ar_cut_0_face(repo, ar_comlpex, i, plane_index));
    }

    // Step 2: handle 1-faces.
    std::vector<std::array<size_t, 3>> subedges;
    subedges.reserve(num_edges);
    for (size_t i = 0; i < num_edges; i++) {
        subedges.push_back(ar_cut_1_face(ar_comlpex, i, plane_index, orientations));
    }

    // Step 3: handle 2-faces.
    std::vector<std::array<size_t, 3>> subfaces;
    subfaces.reserve(num_faces);
    for (size_t i = 0; i < num_faces; i++) {
        subfaces.push_back(ar_cut_2_face(ar_comlpex, i, plane_index, orientations, subedges));
    }

    // Step 4: Remove old faces and update indices.
    {
        std::vector<bool> to_keep(faces.size(), false);
        for (const auto& subface : subfaces) {
            if (subface[0] != INVALID) to_keep[subface[0]] = true;
            if (subface[1] != INVALID) to_keep[subface[1]] = true;
        }

        auto index_map = utils::shrink(faces, [&](size_t i) { return to_keep[i]; });

        // Update face indices in edges.
        for (auto& e : edges) {
            if (e.positive_face != INVALID) e.positive_face = index_map[e.positive_face];
            if (e.negative_face != INVALID) e.negative_face = index_map[e.negative_face];
        }
    }

    // Step 5: consolidate.
    consolidate(ar_comlpex);
    logger().debug("After: {} {} {}",
        ar_comlpex.vertices.size(),
        ar_comlpex.edges.size(),
        ar_comlpex.faces.size());
}

template <typename Scalar>
void add_plane(const PlaneRepo<Scalar, 3>& repo, MIComplex<3>& ar_complex, size_t plane_index)
{
    const size_t num_vertices = ar_complex.vertices.size();
    const size_t num_edges = ar_complex.edges.size();
    const size_t num_faces = ar_complex.faces.size();
    const size_t num_cells = ar_complex.cells.size();
    logger().debug("adding material {}", plane_index);
    logger().debug("Before: {} {} {} {}", num_vertices, num_edges, num_faces, num_cells);

    auto& vertices = ar_complex.vertices;
    auto& edges = ar_complex.edges;
    auto& faces = ar_complex.faces;
    auto& cells = ar_complex.cells;

    // Reserve capacity. TODO: check if this helps.
    vertices.reserve(num_vertices + num_edges);
    edges.reserve(num_edges * 2);
    faces.reserve(num_faces * 2);
    cells.reserve(num_cells * 2);

    // Step 1: handle 0-faces.
    std::vector<int8_t> orientations;
    orientations.reserve(num_vertices);
    for (size_t i = 0; i < num_vertices; i++) {
        orientations.push_back(ar_cut_0_face(repo, ar_complex, i, plane_index));
    }

    // Step 2: handle 1-faces.
    std::vector<std::array<size_t, 3>> subedges;
    subedges.reserve(num_edges);
    for (size_t i = 0; i < num_edges; i++) {
        subedges.push_back(ar_cut_1_face(ar_complex, i, plane_index, orientations));
    }

    // Step 3: handle 2-faces.
    std::vector<std::array<size_t, 3>> subfaces;
    subfaces.reserve(num_faces);
    for (size_t i = 0; i < num_faces; i++) {
        subfaces.push_back(ar_cut_2_face(ar_complex, i, plane_index, orientations, subedges));
    }

    // Step 4: handle 3-faces.
    std::vector<std::array<size_t, 3>> subcells;
    subcells.reserve(num_cells);
    for (size_t i = 0; i < num_cells; i++) {
        subcells.push_back(ar_cut_3_face(ar_complex, i, plane_index, subfaces));
    }

    // Step 5: Remove old cells and update cell indices
    {
        std::vector<bool> to_keep(cells.size(), false);
        for (const auto& subcells : subcells) {
            if (subcell[0] != INVALID) to_keep[subcell[0]] = true;
            if (subcell[1] != INVALID) to_keep[subcell[1]] = true;
        }

        auto index_map = utils::shrink(cells, [&](size_t i) { return to_keep[i]; });

        // Update cell indices in faces.
        for (auto& f : faces) {
            if (f.positive_cell != INVALID) f.positive_cell = index_map[f.positive_cell];
            if (f.negative_cell != INVALID) f.negative_cell = index_map[f.negative_cell];
        }
    }

    // Step 6: consolidate.
    consolidate(ar_complex);
    logger().debug("After: {} {} {} {}",
        ar_complex.vertices.size(),
        ar_complex.edges.size(),
        ar_complex.faces.size(),
        ar_complex.cells.size());
}

} // namespace

namespace simplicial_arrangement::internal {

void add_plane(const PlaneRepo<double, 2>& repo, ARComplex<2>& ar_complex, size_t plane_index)
{
    ::add_plane(repo, ar_complex, plane_index);
}

void add_plane(const PlaneRepo<Int, 2>& repo, ARComplex<2>& ar_complex, size_t plane_index)
{
    ::add_plane(repo, ar_complex, plane_index);
}

void add_plane(const PlaneRepo<double, 3>& repo, ARComplex<3>& ar_complex, size_t plane_index)
{
    ::add_plane(repo, ar_complex, plane_index);
}

void add_plane(const PlaneRepo<Int, 3>& repo, ARComplex<2>& ar_complex, size_t plane_index)
{
    ::add_plane(repo, ar_complex, plane_index);
}

} // namespace simplicial_arrangement::internal
