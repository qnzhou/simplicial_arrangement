#include <simplicial_arrangement/material_interface.h>
#include <spdlog/spdlog.h>
#include <catch2/catch.hpp>

#include <implicit_predicates/implicit_predicates.h>

#include <MaterialRepo.h>
#include <initialize_simplicial_mi_complex.h>

#include "test_mi_cut_utils.h"

namespace {

template <typename Scalar>
void test_2D()
{
    using namespace simplicial_arrangement;
    std::vector<Material<Scalar, 2>> materials;
    materials.push_back({0, 0, 0}); // material 3
    materials.push_back({1, 0, 0}); // material 4
    materials.push_back({0, 1, 0}); // material 5
    materials.push_back({0, 0, 1}); // material 6
    materials.push_back({1, 1, 1}); // material 7

    MaterialRepo<Scalar, 2> repo(materials);
    auto mi_complex = initialize_simplicial_mi_complex<2>();

    SECTION("Tagent case")
    {
        // Inserting material 4.
        size_t material_index = 4;
        auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);
        auto subedges = test_utils::compute_subedges(mi_complex, repo, material_index, orientations);
        auto r = mi_cut_2_face(mi_complex, 0, material_index, orientations, subedges);
        REQUIRE(r[0] == INVALID);
        REQUIRE(r[1] != INVALID);
        REQUIRE(r[2] != INVALID);
    }
    SECTION("No cut")
    {
        size_t material_index = 7;
        auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);
        auto subedges = test_utils::compute_subedges(mi_complex, repo, material_index, orientations);
        auto r = mi_cut_2_face(mi_complex, 0, material_index, orientations, subedges);
        REQUIRE(r[0] == INVALID);
        REQUIRE(r[1] == 0);
        REQUIRE(r[2] == INVALID);

        const auto& f = mi_complex.faces[0];
        REQUIRE(f.material_label == material_index);
    }
    SECTION("Cross cut")
    {
        materials.push_back({-1, 1, 0}); // material 8
        size_t material_index = 8;
        auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);
        auto subedges = test_utils::compute_subedges(mi_complex, repo, material_index, orientations);
        auto r = mi_cut_2_face(mi_complex, 0, material_index, orientations, subedges);
        REQUIRE(r[0] != INVALID);
        REQUIRE(r[1] != INVALID);
        REQUIRE(r[2] != INVALID);
    }
}

template <typename Scalar>
void test_3D()
{
    spdlog::set_level(spdlog::level::info);
    using namespace simplicial_arrangement;
    std::vector<Material<Scalar, 3>> materials;
    materials.push_back({0, 0, 0, 0}); // material 4
    materials.push_back({1, 0, 0, 0}); // material 5
    materials.push_back({0, 1, 0, 0}); // material 6
    materials.push_back({0, 0, 1, 0}); // material 7
    materials.push_back({0, 0, 0, 1}); // material 8
    materials.push_back({1, 1, 1, 1}); // material 9

    MaterialRepo<Scalar, 3> repo(materials);
    auto mi_complex = initialize_simplicial_mi_complex<3>();

    SECTION("Tagent case: point")
    {
        // Inserting material 6.
        size_t material_index = 6;
        auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);
        auto subedges = test_utils::compute_subedges(mi_complex, repo, material_index, orientations);
        auto r = mi_cut_2_face(mi_complex, 0, material_index, orientations, subedges);
        REQUIRE(r[0] == INVALID);
        REQUIRE(r[1] == 0);
        REQUIRE(r[2] != INVALID);
    }
    SECTION("Tagent case: edge")
    {
        materials.push_back({1, 1, 0, 0}); // material 10
        size_t material_index = 10;
        const size_t num_edges = mi_complex.edges.size();
        auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);
        auto subedges = test_utils::compute_subedges(mi_complex, repo, material_index, orientations);
        auto r = mi_cut_2_face(mi_complex, 0, material_index, orientations, subedges);
        REQUIRE(r[0] == INVALID);
        REQUIRE(r[1] == 0);
        REQUIRE(r[2] < num_edges); // Must be one of the existing edges.
    }
    SECTION("Tagent case: face")
    {
        // Inserting material 5.
        size_t material_index = 5;
        auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);
        auto subedges = test_utils::compute_subedges(mi_complex, repo, material_index, orientations);
        auto r = mi_cut_2_face(mi_complex, 0, material_index, orientations, subedges);
        REQUIRE(r[0] == INVALID);
        REQUIRE(r[1] == INVALID);
        REQUIRE(r[2] == INVALID);
    }
    SECTION("No cut")
    {
        size_t material_index = 9;
        auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);
        auto subedges = test_utils::compute_subedges(mi_complex, repo, material_index, orientations);
        auto r = mi_cut_2_face(mi_complex, 0, material_index, orientations, subedges);
        REQUIRE(r[0] == INVALID);
        REQUIRE(r[1] == 0);
        REQUIRE(r[2] == INVALID);
    }
    SECTION("Cross cut")
    {
        materials.push_back({0, -1, 1, 0}); // material 10
        size_t material_index = 10;
        auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);
        auto subedges = test_utils::compute_subedges(mi_complex, repo, material_index, orientations);
        auto r = mi_cut_2_face(mi_complex, 0, material_index, orientations, subedges);
        REQUIRE(r[0] != INVALID);
        REQUIRE(r[1] != INVALID);
        REQUIRE(r[2] != INVALID);

        const auto& f = mi_complex.faces[0];
        const auto& f0 = mi_complex.faces[r[0]];
        const auto& f1 = mi_complex.faces[r[1]];
        REQUIRE(f0.positive_material_label == f.positive_material_label);
        REQUIRE(f0.negative_material_label == f.negative_material_label);
        REQUIRE(f1.positive_material_label == f.positive_material_label);
        REQUIRE(f1.negative_material_label == f.negative_material_label);
    }
}

} // namespace

TEST_CASE("mi_cut_2_face", "[material_interface]")
{
    using namespace simplicial_arrangement;
    SECTION("2D")
    {
#ifndef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
        SECTION("Int") { test_2D<Int>(); }
#endif
        SECTION("double") { test_2D<double>(); }
    }

    SECTION("3D")
    {
#ifndef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
        SECTION("Int") { test_3D<Int>(); }
#endif
        SECTION("double") { test_3D<double>(); }
    }
}
