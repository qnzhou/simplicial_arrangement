#pragma once

#include <ARComplex.h>
#include <PlaneRepo.h>
#include <ar_cut_0_face.h>
#include <ar_cut_1_face.h>
#include <ar_cut_2_face.h>

namespace simplicial_arrangement::test_utils {

template <int DIM>
void initialize_signs(simplicial_arrangement::ARComplex<DIM>& ar_complex, size_t num_materials) {
    if constexpr (DIM == 2) {
        for (auto& f : ar_complex.faces) {
            f.signs.resize(num_materials + DIM + 1, false);
        }
    } else {
        for (auto& c : ar_complex.cells) {
            c.signs.resize(num_materials + DIM + 1, false);
        }
    }
}


template <typename Scalar, int DIM>
std::vector<int8_t> compute_orientations(const simplicial_arrangement::ARComplex<DIM>& ar_complex,
    const simplicial_arrangement::PlaneRepo<Scalar, DIM>& plane_repo,
    size_t plane_index)
{
    const size_t num_vertices = ar_complex.vertices.size();
    std::vector<int8_t> orientations;
    orientations.reserve(num_vertices);
    for (size_t i = 0; i < num_vertices; i++) {
        orientations.push_back(
            simplicial_arrangement::ar_cut_0_face(plane_repo, ar_complex, i, plane_index));
    }
    return orientations;
}

template <typename Scalar, int DIM>
std::vector<std::array<size_t, 3>> compute_subedges(
    simplicial_arrangement::ARComplex<DIM>& ar_complex,
    const simplicial_arrangement::PlaneRepo<Scalar, DIM>& plane_repo,
    size_t plane_index,
    const std::vector<int8_t>& orientations)
{
    const size_t num_edges = ar_complex.edges.size();
    std::vector<std::array<size_t, 3>> subedges;
    subedges.reserve(num_edges);

    for (size_t eid = 0; eid < num_edges; eid++) {
        subedges.push_back(ar_cut_1_face(ar_complex, eid, plane_index, orientations));
    }

    return subedges;
}

template <typename Scalar, int DIM>
std::vector<std::array<size_t, 3>> compute_subfaces(
    simplicial_arrangement::ARComplex<DIM>& ar_complex,
    const simplicial_arrangement::PlaneRepo<Scalar, DIM>& plane_repo,
    size_t plane_index,
    const std::vector<int8_t>& orientations,
    const std::vector<std::array<size_t, 3>>& subedges)
{
    const size_t num_faces = ar_complex.faces.size();
    std::vector<std::array<size_t, 3>> subfaces;
    subfaces.reserve(num_faces);

    for (size_t fid = 0; fid < num_faces; fid++) {
        subfaces.push_back(ar_cut_2_face(ar_complex, fid, plane_index, orientations, subedges));
    }

    return subfaces;
}

} // namespace simplicial_arrangement::test_utils
