#include <simplicial_arrangement/simplicial_arrangement.h>
#include <spdlog/spdlog.h>
#include <catch2/catch.hpp>

#include <PlaneRepo.h>
#include <ar_cut_3_face.h>
#include <initialize_simplicial_ar_complex.h>

#include "test_ar_cut_utils.h"

namespace {

template <typename Scalar>
void test_3D()
{
    spdlog::set_level(spdlog::level::info);
    using namespace simplicial_arrangement;
    std::vector<Plane<Scalar, 3>> planes;
    planes.push_back({0, 0, 0, 0}); // plane 4
    planes.push_back({1, 0, 0, 0}); // plane 5
    planes.push_back({1, -1, -1, -1}); // plane 6
    planes.push_back({1, 1, 1, 0}); // plane 7
    planes.push_back({-1, 0, 0, 1}); // plane 8
    planes.push_back({1, 1, 1, 1}); // plane 9
    planes.push_back({1, 1, -2, -5}); // plane 10
    planes.push_back({-1, 1, -2, 0}); // plane 11

    PlaneRepo<Scalar, 3> repo(planes);
    auto ar_complex = initialize_simplicial_ar_complex<3>();
    test_utils::initialize_signs(ar_complex, planes.size());

    SECTION("Tangent case: cell")
    {
        size_t plane_index = 4;
        auto orientations = test_utils::compute_orientations(ar_complex, repo, plane_index);
        auto subedges =
            test_utils::compute_subedges(ar_complex, repo, plane_index, orientations);
        auto subfaces =
            test_utils::compute_subfaces(ar_complex, repo, plane_index, orientations, subedges);
        auto r = ar_cut_3_face(ar_complex, 0, plane_index, subfaces);
        REQUIRE(r[0] == INVALID);
        REQUIRE(r[1] == INVALID);
        REQUIRE(r[2] == INVALID);
    }
    SECTION("Tangent case: face")
    {
        size_t plane_index = 5;
        auto orientations = test_utils::compute_orientations(ar_complex, repo, plane_index);
        auto subedges =
            test_utils::compute_subedges(ar_complex, repo, plane_index, orientations);
        auto subfaces =
            test_utils::compute_subfaces(ar_complex, repo, plane_index, orientations, subedges);
        auto r = ar_cut_3_face(ar_complex, 0, plane_index, subfaces);
        REQUIRE(r[0] == 0);
        REQUIRE(r[1] == INVALID);
        REQUIRE(r[2] != INVALID);

        const auto& c0 = ar_complex.cells[r[0]];
        REQUIRE(c0.signs[plane_index]);
    }
    SECTION("Tangent case: point")
    {
        size_t plane_index = 7;
        auto orientations = test_utils::compute_orientations(ar_complex, repo, plane_index);
        auto subedges =
            test_utils::compute_subedges(ar_complex, repo, plane_index, orientations);
        auto subfaces =
            test_utils::compute_subfaces(ar_complex, repo, plane_index, orientations, subedges);
        auto r = ar_cut_3_face(ar_complex, 0, plane_index, subfaces);
        REQUIRE(r[0] == 0);
        REQUIRE(r[1] == INVALID);
        REQUIRE(r[2] == INVALID);

        const auto& c0 = ar_complex.cells[r[0]];
        REQUIRE(c0.signs[plane_index]);
    }
    SECTION("No cut")
    {
        size_t plane_index = 9;
        auto orientations = test_utils::compute_orientations(ar_complex, repo, plane_index);
        auto subedges =
            test_utils::compute_subedges(ar_complex, repo, plane_index, orientations);
        auto subfaces =
            test_utils::compute_subfaces(ar_complex, repo, plane_index, orientations, subedges);
        auto r = ar_cut_3_face(ar_complex, 0, plane_index, subfaces);
        REQUIRE(r[0] == 0);
        REQUIRE(r[1] == INVALID);
        REQUIRE(r[2] == INVALID);

        const auto& c0 = ar_complex.cells[r[0]];
        REQUIRE(c0.signs[plane_index]);
    }
    SECTION("Cross cut")
    {
        SECTION("Case 1: cut through edge") {
            size_t plane_index = 8;
            auto orientations = test_utils::compute_orientations(ar_complex, repo, plane_index);
            auto subedges =
                test_utils::compute_subedges(ar_complex, repo, plane_index, orientations);
            auto subfaces =
                test_utils::compute_subfaces(ar_complex, repo, plane_index, orientations, subedges);
            auto r = ar_cut_3_face(ar_complex, 0, plane_index, subfaces);
            REQUIRE(r[0] != INVALID);
            REQUIRE(r[1] != INVALID);
            REQUIRE(r[2] != INVALID);

            const auto& c0 = ar_complex.cells[r[0]];
            const auto& c1 = ar_complex.cells[r[1]];
            REQUIRE(c0.signs[plane_index]);
            REQUIRE(!c1.signs[plane_index]);

            const auto& f = ar_complex.faces[r[2]];
            REQUIRE(f.supporting_plane == plane_index);
            REQUIRE(f.edges.size() == 3);
        }
        SECTION("Case 2: triangle cross section") {
            size_t plane_index = 6;
            auto orientations = test_utils::compute_orientations(ar_complex, repo, plane_index);
            auto subedges =
                test_utils::compute_subedges(ar_complex, repo, plane_index, orientations);
            auto subfaces =
                test_utils::compute_subfaces(ar_complex, repo, plane_index, orientations, subedges);
            auto r = ar_cut_3_face(ar_complex, 0, plane_index, subfaces);
            REQUIRE(r[0] != INVALID);
            REQUIRE(r[1] != INVALID);
            REQUIRE(r[2] != INVALID);

            const auto& c0 = ar_complex.cells[r[0]];
            const auto& c1 = ar_complex.cells[r[1]];
            REQUIRE(c0.signs[plane_index]);
            REQUIRE(!c1.signs[plane_index]);

            const auto& f = ar_complex.faces[r[2]];
            REQUIRE(f.supporting_plane == plane_index);
            REQUIRE(f.edges.size() == 3);
        }
        SECTION("Case 3: quad cross section") {
            size_t plane_index = 10;
            auto orientations = test_utils::compute_orientations(ar_complex, repo, plane_index);
            auto subedges =
                test_utils::compute_subedges(ar_complex, repo, plane_index, orientations);
            auto subfaces =
                test_utils::compute_subfaces(ar_complex, repo, plane_index, orientations, subedges);
            auto r = ar_cut_3_face(ar_complex, 0, plane_index, subfaces);
            REQUIRE(r[0] != INVALID);
            REQUIRE(r[1] != INVALID);
            REQUIRE(r[2] != INVALID);

            const auto& c0 = ar_complex.cells[r[0]];
            const auto& c1 = ar_complex.cells[r[1]];
            REQUIRE(c0.signs[plane_index]);
            REQUIRE(!c1.signs[plane_index]);

            const auto& f = ar_complex.faces[r[2]];
            REQUIRE(f.supporting_plane == plane_index);
            REQUIRE(f.edges.size() == 4);
        }
        SECTION("Case 4: cut through vertex") {
            size_t plane_index = 11;
            auto orientations = test_utils::compute_orientations(ar_complex, repo, plane_index);
            auto subedges =
                test_utils::compute_subedges(ar_complex, repo, plane_index, orientations);
            auto subfaces =
                test_utils::compute_subfaces(ar_complex, repo, plane_index, orientations, subedges);
            auto r = ar_cut_3_face(ar_complex, 0, plane_index, subfaces);
            REQUIRE(r[0] != INVALID);
            REQUIRE(r[1] != INVALID);
            REQUIRE(r[2] != INVALID);

            const auto& c0 = ar_complex.cells[r[0]];
            const auto& c1 = ar_complex.cells[r[1]];
            REQUIRE(c0.signs[plane_index]);
            REQUIRE(!c1.signs[plane_index]);

            const auto& f = ar_complex.faces[r[2]];
            REQUIRE(f.supporting_plane == plane_index);
            REQUIRE(f.edges.size() == 3);
        }
    }
}

} // namespace

TEST_CASE("ar_cut_3_face", "[plane_interface]")
{
    using namespace simplicial_arrangement;

    SECTION("3D")
    {
        SECTION("Int") { test_3D<Int>(); }
        SECTION("double") { test_3D<double>(); }
    }
}
