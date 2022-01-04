#include <simplicial_arrangement/simplicial_arrangement.h>
#include <spdlog/spdlog.h>
#include <catch2/catch.hpp>

#include <implicit_predicates/implicit_predicates.h>

#include <PlaneRepo.h>
#include <initialize_simplicial_ar_complex.h>

#include "test_ar_cut_utils.h"

namespace {

template <typename Scalar>
void test_2D()
{
    using namespace simplicial_arrangement;
    std::vector<Plane<Scalar, 2>> planes;
    planes.push_back({0, 0, 0}); // plane 3
    planes.push_back({1, 0, 0}); // plane 4
    planes.push_back({0, 1, 0}); // plane 5
    planes.push_back({0, 0, 1}); // plane 6
    planes.push_back({1, 1, 1}); // plane 7
    planes.push_back({-1, 1, 1}); // plane 8

    PlaneRepo<Scalar, 2> repo(planes);
    auto ar_complex = initialize_simplicial_ar_complex<2>(planes.size() * 2);

    SECTION("Tangent case")
    {
        size_t plane_index = 4;
        auto orientations = test_utils::compute_orientations(ar_complex, repo, plane_index);
        REQUIRE(ar_complex.vertices.size() == 3);
        REQUIRE(ar_complex.edges.size() == 3);
        REQUIRE(ar_complex.faces.size() == 1);

        const size_t num_edges = ar_complex.edges.size();
        for (size_t eid = 0; eid < num_edges; eid++) {
            const auto e = ar_complex.edges[eid];
            auto r = ar_cut_1_face(ar_complex, eid, plane_index, orientations);
            if (e.vertices[0] == 0 || e.vertices[1] == 0) {
                REQUIRE(r[0] == eid);
                REQUIRE(r[1] == INVALID);
                REQUIRE(r[2] != INVALID);
            } else {
                REQUIRE(r[0] == INVALID);
                REQUIRE(r[1] == INVALID);
            }
        }

        REQUIRE(ar_complex.vertices.size() == 3);
        REQUIRE(ar_complex.edges.size() == 3);
        REQUIRE(ar_complex.faces.size() == 1);
    }

    SECTION("Non tangent case")
    {
        size_t plane_index = 8;
        auto orientations = test_utils::compute_orientations(ar_complex, repo, plane_index);
        REQUIRE(ar_complex.vertices.size() == 3);
        REQUIRE(ar_complex.edges.size() == 3);

        const size_t num_edges = ar_complex.edges.size();
        for (size_t eid = 0; eid < num_edges; eid++) {
            const auto e = ar_complex.edges[eid];
            auto r = ar_cut_1_face(ar_complex, eid, plane_index, orientations);
            if (e.vertices[0] == 0 || e.vertices[1] == 0) {
                // Cross cut.
                REQUIRE(r[0] != INVALID);
                REQUIRE(r[1] != INVALID);
                REQUIRE(r[2] != INVALID);

                const auto& e0 = ar_complex.edges[r[0]];
                const auto& e1 = ar_complex.edges[r[1]];
                REQUIRE(e0.supporting_plane == e.supporting_plane);
                REQUIRE(e1.supporting_plane == e.supporting_plane);
            } else {
                // Miss.
                REQUIRE(r[0] == eid);
                REQUIRE(r[1] == INVALID);
                REQUIRE(r[2] == INVALID);
            }
        }

        REQUIRE(ar_complex.vertices.size() == 5); // 3 old vertices, 2 new vertex.
        REQUIRE(ar_complex.edges.size() == 7); // 3 old edges, 4 new edges.
    }
}

template <typename Scalar>
void test_3D()
{
    using namespace simplicial_arrangement;
    std::vector<Plane<Scalar, 3>> planes;
    planes.push_back({0, 0, 0, 0}); // plane 4
    planes.push_back({1, 0, 0, 0}); // plane 5
    planes.push_back({1, 0, 0, -1}); // plane 6
    planes.push_back({0, 0, 1, 0}); // plane 7
    planes.push_back({0, 0, 0, 1}); // plane 8
    planes.push_back({2, 2, 2, 2}); // plane 9
    planes.push_back({1, -1, -1, -1}); // plane 10

    PlaneRepo<Scalar, 3> repo(planes);
    auto ar_complex = initialize_simplicial_ar_complex<3>(planes.size() * 2);

    SECTION("Tangent case") {
        size_t plane_index = 5;
        auto orientations = test_utils::compute_orientations(ar_complex, repo, plane_index);

        const size_t num_edges = ar_complex.edges.size();
        for (size_t eid = 0; eid < num_edges; eid++) {
            const auto e = ar_complex.edges[eid];
            const auto r = ar_cut_1_face(ar_complex, eid, plane_index, orientations);
            if (e.vertices[0] == 0 || e.vertices[1] == 0) {
                // Tangent at vertex.
                REQUIRE(r[0] != INVALID);
                REQUIRE(r[1] == INVALID);
                REQUIRE(r[2] != INVALID);
                REQUIRE(r[2] != 0);
            } else {
                // Edge and plane are coplanar.
                REQUIRE(r[0] == INVALID);
                REQUIRE(r[1] == INVALID);
            }
        }
    }

    SECTION("Tangent case 2") {
        size_t plane_index = 6;
        auto orientations = test_utils::compute_orientations(ar_complex, repo, plane_index);

        SECTION("edge 0")
        {
            const auto r = ar_cut_1_face(ar_complex, 0, plane_index, orientations);
            REQUIRE(r[0] == 0);
            REQUIRE(r[1] == INVALID);
            REQUIRE(r[2] == 1);
        }
        SECTION("edge 1")
        {
            const auto r = ar_cut_1_face(ar_complex, 1, plane_index, orientations);
            REQUIRE(r[0] == 1);
            REQUIRE(r[1] == INVALID);
            REQUIRE(r[2] == 2);
        }
        SECTION("edge 2")
        {
            const auto r = ar_cut_1_face(ar_complex, 2, plane_index, orientations);
            REQUIRE(r[0] != INVALID);
            REQUIRE(r[1] != INVALID);
            REQUIRE(r[2] == 4);
        }
        SECTION("edge 3")
        {
            const auto r = ar_cut_1_face(ar_complex, 3, plane_index, orientations);
            REQUIRE(r[0] == INVALID);
            REQUIRE(r[1] == INVALID);
        }
        SECTION("edge 4")
        {
            const auto r = ar_cut_1_face(ar_complex, 4, plane_index, orientations);
            REQUIRE(r[0] == INVALID);
            REQUIRE(r[1] == 4);
            REQUIRE(r[2] == 1);
        }
        SECTION("edge 5")
        {
            const auto r = ar_cut_1_face(ar_complex, 5, plane_index, orientations);
            REQUIRE(r[0] == INVALID);
            REQUIRE(r[1] == 5);
            REQUIRE(r[2] == 2);
        }
    }

    SECTION("Non-tangent case") {
        size_t plane_index = 10;
        auto orientations = test_utils::compute_orientations(ar_complex, repo, plane_index);
        REQUIRE(ar_complex.vertices.size() == 4);
        REQUIRE(ar_complex.edges.size() == 6);

        const size_t num_edges = ar_complex.edges.size();
        for (size_t eid = 0; eid < num_edges; eid++) {
            const auto e = ar_complex.edges[eid];
            const auto r = ar_cut_1_face(ar_complex, eid, plane_index, orientations);
            if (e.vertices[0] == 0 || e.vertices[1] == 0) {
                // Cross cut.
                REQUIRE(r[0] != INVALID);
                REQUIRE(r[1] != INVALID);
                REQUIRE(r[2] != INVALID);

                const auto& e0 = ar_complex.edges[r[0]];
                const auto& e1 = ar_complex.edges[r[1]];
                REQUIRE(e.supporting_planes == e0.supporting_planes);
                REQUIRE(e.supporting_planes == e1.supporting_planes);
            } else {
                // Miss.
                REQUIRE(r[0] == INVALID);
                REQUIRE(r[1] == eid);
                REQUIRE(r[2] == INVALID);
            }
        }
        REQUIRE(ar_complex.vertices.size() == 7); // 4 old vertices, 3 new vertices.
        REQUIRE(ar_complex.edges.size() == 12); // 6 old edges, 6 new edges.
    }
}
}

TEST_CASE("ar_cut_1_face", "[arrangement]")
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
