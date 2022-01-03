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

    SECTION("Tangent case")
    {
        size_t material_index = 5;
        auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);

        const size_t num_edges = mi_complex.edges.size();
        for (size_t eid = 0; eid < num_edges; eid++) {
            const auto e = mi_complex.edges[eid];
            if (e.vertices[0] == 0 && e.vertices[1] == 1) {
                // Cut at vertex.
                auto r = mi_cut_1_face(mi_complex, eid, material_index, orientations);
                REQUIRE(r[2] == 0);
                REQUIRE(r[0] == INVALID);
                REQUIRE(r[1] == eid);
            } else if (e.vertices[0] == 2 && e.vertices[1] == 0) {
                // Cut at edge.
                auto r = mi_cut_1_face(mi_complex, eid, material_index, orientations);
                REQUIRE(r[2] != INVALID);
                REQUIRE(r[0] == INVALID);
                REQUIRE(r[1] == INVALID);
            }
        }
    }

    SECTION("Non-tangent case")
    {
        // Add edge mid point.
        mi_complex.vertices.push_back({0, 5, 6}); // Vertex 3

        MIEdge<2> cut_edge;
        cut_edge.vertices = {0, 3};
        cut_edge.positive_material_label = 5;
        cut_edge.negative_material_label = 6;
        mi_complex.edges.push_back(std::move(cut_edge));
        const auto e = mi_complex.edges.back();

        SECTION("Cross cut")
        {
            size_t material_index = 4;
            auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);
            auto r = mi_cut_1_face(
                mi_complex, mi_complex.edges.size() - 1, material_index, orientations);
            REQUIRE(r[0] != INVALID);
            REQUIRE(r[1] != INVALID);
            REQUIRE(r[2] != INVALID);

            const auto& e0 = mi_complex.edges[r[0]];
            const auto& e1 = mi_complex.edges[r[1]];
            REQUIRE(e0.vertices[1] == 3);
            REQUIRE(e0.vertices[0] == e1.vertices[1]);
            REQUIRE(e1.vertices[0] == 0);

            REQUIRE(e.positive_material_label == e0.positive_material_label);
            REQUIRE(e.negative_material_label == e0.negative_material_label);
            REQUIRE(e.positive_material_label == e1.positive_material_label);
            REQUIRE(e.negative_material_label == e1.negative_material_label);

            auto& v = mi_complex.vertices[r[2]];
            std::sort(v.begin(), v.end());
            REQUIRE(v == Joint<2>({4, 5, 6}));
        }

        SECTION("No cut")
        {
            size_t material_index = 7;
            auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);
            auto r = mi_cut_1_face(
                mi_complex, mi_complex.edges.size() - 1, material_index, orientations);
            REQUIRE(r[0] == INVALID);
            REQUIRE(r[1] == mi_complex.edges.size() - 1);
            REQUIRE(r[2] == INVALID);
        }
    }
}

template <typename Scalar>
void test_3D()
{
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

    SECTION("Tangent case")
    {
        size_t material_index = 5;
        auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);

        const size_t num_edges = mi_complex.edges.size();
        for (size_t eid = 0; eid < num_edges; eid++) {
            const auto e = mi_complex.edges[eid];
            if (e.vertices[0] == 0 && e.vertices[1] == 1) {
                // Cut at vertex.
                auto r = mi_cut_1_face(mi_complex, eid, material_index, orientations);
                REQUIRE(r[2] == 1);
                REQUIRE(r[0] == INVALID);
                REQUIRE(r[1] == eid);
            } else if (e.vertices[0] == 2 && e.vertices[1] == 3) {
                // Cut at edge.
                auto r = mi_cut_1_face(mi_complex, eid, material_index, orientations);
                REQUIRE(r[2] != INVALID);
                REQUIRE(r[0] == INVALID);
                REQUIRE(r[1] == INVALID);
            }
        }
    }

    SECTION("Non-tangent case")
    {
        mi_complex.vertices.push_back({3, 5, 6, 7}); // vertex 4
        MIEdge<3> cut_edge;
        cut_edge.vertices = {3, 4};
        cut_edge.supporting_materials = {5, 6, 7};
        mi_complex.edges.push_back(std::move(cut_edge)); // edge 6;

        SECTION("Cross cut")
        {
            const auto e = mi_complex.edges.back();
            size_t material_index = 8;
            auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);
            auto r = mi_cut_1_face(mi_complex, 6, material_index, orientations);
            REQUIRE(r[0] != INVALID);
            REQUIRE(r[1] != INVALID);
            REQUIRE(r[2] != INVALID);

            const auto& e0 = mi_complex.edges[r[0]];
            const auto& e1 = mi_complex.edges[r[1]];
            REQUIRE(e0.vertices[1] == 4);
            REQUIRE(e0.vertices[0] == e1.vertices[1]);
            REQUIRE(e1.vertices[0] == 3);

            REQUIRE(e.supporting_materials == e0.supporting_materials);
            REQUIRE(e.supporting_materials == e1.supporting_materials);

            auto& v = mi_complex.vertices[r[2]];
            std::sort(v.begin(), v.end());
            REQUIRE(v == Joint<3>({5, 6, 7, 8}));
        }

        SECTION("No cut")
        {
            size_t material_index = 9;
            auto orientations = test_utils::compute_orientations(mi_complex, repo, material_index);
            auto r = mi_cut_1_face(mi_complex, 6, material_index, orientations);
            REQUIRE(r[0] == INVALID);
            REQUIRE(r[1] == 6);
            REQUIRE(r[2] == INVALID);
        }
    }
}

} // namespace

TEST_CASE("mi_cut_1_face", "[material_interface]")
{
    using namespace simplicial_arrangement;
    SECTION("2D")
    {
        SECTION("Int") { test_2D<Int>(); }
        SECTION("double") { test_2D<double>(); }
    }

    SECTION("3D")
    {
        SECTION("Int") { test_3D<Int>(); }
        SECTION("double") { test_3D<double>(); }
    }
}
