#pragma once

#include <MIComplex.h>
#include <MaterialRepo.h>
#include <mi_cut_0_face.h>
#include <mi_cut_1_face.h>
#include <mi_cut_2_face.h>

namespace simplicial_arrangement::test_utils {

template <typename Scalar, int DIM>
std::vector<int8_t> compute_orientations(const simplicial_arrangement::MIComplex<DIM>& mi_complex,
    const simplicial_arrangement::MaterialRepo<Scalar, DIM>& material_repo,
    size_t material_index)
{
    const size_t num_vertices = mi_complex.vertices.size();
    std::vector<int8_t> orientations;
    orientations.reserve(num_vertices);
    for (size_t i = 0; i < num_vertices; i++) {
        orientations.push_back(
            simplicial_arrangement::mi_cut_0_face(material_repo, mi_complex, i, material_index));
    }
    return orientations;
}

template <typename Scalar, int DIM>
std::vector<std::array<size_t, 3>> compute_subedges(
    simplicial_arrangement::MIComplex<DIM>& mi_complex,
    const simplicial_arrangement::MaterialRepo<Scalar, DIM>& material_repo,
    size_t material_index,
    const std::vector<int8_t>& orientations)
{
    const size_t num_edges = mi_complex.edges.size();
    std::vector<std::array<size_t, 3>> subedges;
    subedges.reserve(num_edges);

    for (size_t eid = 0; eid < num_edges; eid++) {
        subedges.push_back(mi_cut_1_face(mi_complex, eid, material_index, orientations));
    }

    return subedges;
}

template <typename Scalar, int DIM>
std::vector<std::array<size_t, 3>> compute_subfaces(
    simplicial_arrangement::MIComplex<DIM>& mi_complex,
    const simplicial_arrangement::MaterialRepo<Scalar, DIM>& material_repo,
    size_t material_index,
    const std::vector<int8_t>& orientations,
    const std::vector<std::array<size_t, 3>>& subedges)
{
    const size_t num_faces = mi_complex.faces.size();
    std::vector<std::array<size_t, 3>> subfaces;
    subfaces.reserve(num_faces);

    for (size_t fid = 0; fid < num_faces; fid++) {
        subfaces.push_back(mi_cut_2_face(mi_complex, fid, material_index, orientations, subedges));
    }

    return subfaces;
}

} // namespace simplicial_arrangement::test_utils
