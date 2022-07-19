#include <simplicial_arrangement/material_interface.h>
#include <spdlog/spdlog.h>
#include <catch2/catch.hpp>

#include <implicit_predicates/implicit_predicates.h>

#include <MaterialRepo.h>
#include <initialize_simplicial_mi_complex.h>
#include <mi_cut_3_face.h>

#include "test_mi_cut_utils.h"

namespace {

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

    SECTION("Tangent case: cell")
    {
        materials.push_back({0, 0, 0, 0}); // material 10
        size_t material_index = 10;
        auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);
        auto subedges =
            test_utils::compute_subedges(mi_complex, repo, material_index, orientations);
        auto subfaces =
            test_utils::compute_subfaces(mi_complex, repo, material_index, orientations, subedges);
        auto r = mi_cut_3_face(mi_complex, 0, material_index, subfaces);
        REQUIRE(r[0] == INVALID);
        REQUIRE(r[1] == INVALID);
        REQUIRE(r[2] == INVALID);
    }
    SECTION("Tangent case: face")
    {
        size_t material_index = 5;
        auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);
        auto subedges =
            test_utils::compute_subedges(mi_complex, repo, material_index, orientations);
        auto subfaces =
            test_utils::compute_subfaces(mi_complex, repo, material_index, orientations, subedges);
        auto r = mi_cut_3_face(mi_complex, 0, material_index, subfaces);
        REQUIRE(r[0] == INVALID);
        REQUIRE(r[1] == 0);
        REQUIRE(r[2] != INVALID);
    }
    SECTION("Tangent case: edge")
    {
        materials.push_back({1, 1, 0, 0}); // material 10
        size_t material_index = 10;
        auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);
        auto subedges =
            test_utils::compute_subedges(mi_complex, repo, material_index, orientations);
        auto subfaces =
            test_utils::compute_subfaces(mi_complex, repo, material_index, orientations, subedges);
        auto r = mi_cut_3_face(mi_complex, 0, material_index, subfaces);
        REQUIRE(r[0] == INVALID);
        REQUIRE(r[1] == 0);
        REQUIRE(r[2] == INVALID);
    }
    SECTION("Tangent case: point")
    {
        materials.push_back({1, 1, 1, 0}); // material 10
        size_t material_index = 10;
        auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);
        auto subedges =
            test_utils::compute_subedges(mi_complex, repo, material_index, orientations);
        auto subfaces =
            test_utils::compute_subfaces(mi_complex, repo, material_index, orientations, subedges);
        auto r = mi_cut_3_face(mi_complex, 0, material_index, subfaces);
        REQUIRE(r[0] == INVALID);
        REQUIRE(r[1] == 0);
        REQUIRE(r[2] == INVALID);
    }
    SECTION("No cut")
    {
        size_t material_index = 9;
        auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);
        auto subedges =
            test_utils::compute_subedges(mi_complex, repo, material_index, orientations);
        auto subfaces =
            test_utils::compute_subfaces(mi_complex, repo, material_index, orientations, subedges);
        auto r = mi_cut_3_face(mi_complex, 0, material_index, subfaces);
        REQUIRE(r[0] == INVALID);
        REQUIRE(r[1] == 0);
        REQUIRE(r[2] == INVALID);
    }
    SECTION("Cross cut")
    {
        SECTION("Case 1: cut through edge")
        {
            materials.push_back({1, -1, 0, 0}); // material 10
            size_t material_index = 10;
            auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);
            auto subedges =
                test_utils::compute_subedges(mi_complex, repo, material_index, orientations);
            auto subfaces = test_utils::compute_subfaces(
                mi_complex, repo, material_index, orientations, subedges);
            auto r = mi_cut_3_face(mi_complex, 0, material_index, subfaces);
            REQUIRE(r[0] != INVALID);
            REQUIRE(r[1] != INVALID);
            REQUIRE(r[2] != INVALID);

            const auto& c = mi_complex.cells[0];
            const auto& c0 = mi_complex.cells[r[0]];
            const auto& c1 = mi_complex.cells[r[1]];
            REQUIRE(c0.material_label == c.material_label);
            REQUIRE(c1.material_label == material_index);

            const auto& f = mi_complex.faces[r[2]];
            REQUIRE(f.edges.size() == 3);
        }
        SECTION("Case 2: triangle cross section")
        {
            materials.push_back({1, -1, -1, -1}); // material 10
            size_t material_index = 10;
            auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);
            auto subedges =
                test_utils::compute_subedges(mi_complex, repo, material_index, orientations);
            auto subfaces = test_utils::compute_subfaces(
                mi_complex, repo, material_index, orientations, subedges);
            auto r = mi_cut_3_face(mi_complex, 0, material_index, subfaces);
            REQUIRE(r[0] != INVALID);
            REQUIRE(r[1] != INVALID);
            REQUIRE(r[2] != INVALID);

            const auto& c = mi_complex.cells[0];
            const auto& c0 = mi_complex.cells[r[0]];
            const auto& c1 = mi_complex.cells[r[1]];
            REQUIRE(c0.material_label == c.material_label);
            REQUIRE(c1.material_label == material_index);

            const auto& f = mi_complex.faces[r[2]];
            REQUIRE(f.edges.size() == 3);
        }
        SECTION("Case 3: qaud cross section")
        {
            materials.push_back({1, 1, -2, -3}); // material 10
            size_t material_index = 10;
            auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);
            auto subedges =
                test_utils::compute_subedges(mi_complex, repo, material_index, orientations);
            auto subfaces = test_utils::compute_subfaces(
                mi_complex, repo, material_index, orientations, subedges);
            auto r = mi_cut_3_face(mi_complex, 0, material_index, subfaces);
            REQUIRE(r[0] != INVALID);
            REQUIRE(r[1] != INVALID);
            REQUIRE(r[2] != INVALID);

            const auto& c = mi_complex.cells[0];
            const auto& c0 = mi_complex.cells[r[0]];
            const auto& c1 = mi_complex.cells[r[1]];
            REQUIRE(c0.material_label == c.material_label);
            REQUIRE(c1.material_label == material_index);

            const auto& f = mi_complex.faces[r[2]];
            REQUIRE(f.edges.size() == 4);
        }
        SECTION("Case 4: cut through vertex")
        {
            materials.push_back({1, 1, -2, 0}); // material 10
            size_t material_index = 10;
            auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);
            auto subedges =
                test_utils::compute_subedges(mi_complex, repo, material_index, orientations);
            auto subfaces = test_utils::compute_subfaces(
                mi_complex, repo, material_index, orientations, subedges);
            auto r = mi_cut_3_face(mi_complex, 0, material_index, subfaces);
            REQUIRE(r[0] != INVALID);
            REQUIRE(r[1] != INVALID);
            REQUIRE(r[2] != INVALID);

            const auto& c = mi_complex.cells[0];
            const auto& c0 = mi_complex.cells[r[0]];
            const auto& c1 = mi_complex.cells[r[1]];
            REQUIRE(c0.material_label == c.material_label);
            REQUIRE(c1.material_label == material_index);

            const auto& f = mi_complex.faces[r[2]];
            REQUIRE(f.edges.size() == 3);
        }
    }
}

} // namespace

TEST_CASE("mi_cut_3_face", "[material_interface]")
{
    using namespace simplicial_arrangement;

    SECTION("3D")
    {
#ifndef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
        SECTION("Int") { test_3D<Int>(); }
#endif
        SECTION("double") { test_3D<double>(); }
    }
}
