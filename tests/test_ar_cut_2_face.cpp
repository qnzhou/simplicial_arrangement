#include <simplicial_arrangement/simplicial_arrangement.h>
#include <spdlog/spdlog.h>
#include <catch2/catch.hpp>

#include <PlaneRepo.h>
#include <initialize_simplicial_ar_complex.h>

#include "test_ar_cut_utils.h"

namespace {

template <int DIM>
void check_face_orientation(const simplicial_arrangement::ARComplex<DIM>& ar_complex,
    size_t positive_fid,
    size_t negative_fid,
    size_t cut_edge_eid)
{
    using namespace simplicial_arrangement;
    assert(positive_fid != INVALID);
    assert(negative_fid != INVALID);
    assert(cut_edge_eid != INVALID);

    const auto& f0 = ar_complex.faces[positive_fid];
    const auto& f1 = ar_complex.faces[negative_fid];

    auto get_shared_vertex = [](const auto& e0, const auto& e1) {
        if (e0.vertices[0] == e1.vertices[0] || e0.vertices[0] == e1.vertices[1]) {
            return e0.vertices[0];
        } else {
            return e0.vertices[1];
        }
    };

    auto check_orientation = [&](const auto& f) {
        bool r = false;
        for (size_t i = 0; i < f.edges.size(); i++) {
            const size_t curr_eid = f.edges[i];
            const size_t next_eid = f.edges[(i + 1) % f.edges.size()];
            const auto& curr_e = ar_complex.edges[curr_eid];
            const auto& next_e = ar_complex.edges[next_eid];
            REQUIRE((curr_e.vertices[0] == next_e.vertices[0] ||
                     curr_e.vertices[0] == next_e.vertices[1] ||
                     curr_e.vertices[1] == next_e.vertices[0] ||
                     curr_e.vertices[1] == next_e.vertices[1]));

            if (curr_eid == cut_edge_eid) {
                if (get_shared_vertex(curr_e, next_e) == curr_e.vertices[1]) {
                    r = true;
                } else {
                    r = false;
                }
            }
        }
        return r;
    };

    bool cut_edge_consistent_with_positive_subface = check_orientation(f0);
    bool cut_edge_consistent_with_negative_subface = check_orientation(f1);
    REQUIRE(cut_edge_consistent_with_positive_subface);
    REQUIRE(!cut_edge_consistent_with_negative_subface);
}

template <typename Scalar>
void test_2D()
{
    using namespace simplicial_arrangement;
    std::vector<Plane<Scalar, 2>> planes;
    planes.push_back({0, 0, 0}); // plane 3
    planes.push_back({1, 0, 0}); // plane 4
    planes.push_back({0, 1, 0}); // plane 5
    planes.push_back({0, 0, 1}); // plane 6
    planes.push_back({1, 2, 3}); // plane 7
    planes.push_back({1, -1, -1}); // plane 8

    PlaneRepo<Scalar, 2> repo(planes);
    auto ar_complex = initialize_simplicial_ar_complex<2>(planes.size() * 2);
    test_utils::initialize_signs(ar_complex, planes.size());

    SECTION("Tagent case")
    {
        size_t plane_index = 4;
        auto orientations = test_utils::compute_orientations(ar_complex, repo, plane_index);
        auto subedges = test_utils::compute_subedges(ar_complex, repo, plane_index, orientations);
        auto r = ar_cut_2_face(ar_complex, 0, plane_index, orientations, subedges);
        REQUIRE(r[0] == 0);
        REQUIRE(r[1] == INVALID);
        REQUIRE(r[2] != INVALID);
        REQUIRE(ar_complex.faces[0].signs[plane_index]);
    }

    SECTION("No cut")
    {
        size_t plane_index = 7;
        auto orientations = test_utils::compute_orientations(ar_complex, repo, plane_index);
        auto subedges = test_utils::compute_subedges(ar_complex, repo, plane_index, orientations);
        auto r = ar_cut_2_face(ar_complex, 0, plane_index, orientations, subedges);
        REQUIRE(r[0] == 0);
        REQUIRE(r[1] == INVALID);
        REQUIRE(r[2] == INVALID);
        REQUIRE(ar_complex.faces[0].signs[plane_index]);
    }

    SECTION("Cross cut")
    {
        size_t plane_index = 8;
        auto orientations = test_utils::compute_orientations(ar_complex, repo, plane_index);
        auto subedges = test_utils::compute_subedges(ar_complex, repo, plane_index, orientations);
        auto r = ar_cut_2_face(ar_complex, 0, plane_index, orientations, subedges);
        REQUIRE(r[0] != INVALID);
        REQUIRE(r[1] != INVALID);
        REQUIRE(r[2] != INVALID);

        REQUIRE(ar_complex.vertices.size() == 5); // 3 old + 2 new.
        REQUIRE(ar_complex.edges.size() == 8); // 3 old + 5 new.
        REQUIRE(ar_complex.faces.size() == 3); // 1 old + 2 new.

        const auto& f0 = ar_complex.faces[r[0]];
        const auto& f1 = ar_complex.faces[r[1]];

        REQUIRE(f0.signs[plane_index]);
        REQUIRE(!f1.signs[plane_index]);
        REQUIRE(f0.edges.size() == 3);
        REQUIRE(f1.edges.size() == 4);

        const auto& e = ar_complex.edges[r[2]];
        REQUIRE(e.supporting_plane == plane_index);

        check_face_orientation(ar_complex, r[0], r[1], r[2]);
    }
}

template <typename Scalar>
void test_3D()
{
    spdlog::set_level(spdlog::level::info);
    using namespace simplicial_arrangement;
    std::vector<Plane<Scalar, 3>> planes;
    planes.push_back({-1, -1, -1, 0}); // plane 4
    planes.push_back({1, 0, 0, 0}); // plane 5
    planes.push_back({0, 1, 0, 0}); // plane 6
    planes.push_back({0, 0, 1, 0}); // plane 7
    planes.push_back({0, -1, -1, 1}); // plane 8
    planes.push_back({1, 1, 1, 1}); // plane 9

    PlaneRepo<Scalar, 3> repo(planes);
    auto ar_complex = initialize_simplicial_ar_complex<3>(planes.size() * 2);
    test_utils::initialize_signs(ar_complex, planes.size());

    SECTION("Tagent case: vertex")
    {
        // Inserting plane 4.
        size_t plane_index = 4;
        auto orientations = test_utils::compute_orientations(ar_complex, repo, plane_index);
        auto subedges = test_utils::compute_subedges(ar_complex, repo, plane_index, orientations);
        auto r = ar_cut_2_face(ar_complex, 0, plane_index, orientations, subedges);
        REQUIRE(r[0] == INVALID);
        REQUIRE(r[1] == 0);
        REQUIRE(r[2] == INVALID);
    }
    SECTION("Tagent case: edge")
    {
        // Inserting plane 6.
        size_t plane_index = 6;
        auto orientations = test_utils::compute_orientations(ar_complex, repo, plane_index);
        auto subedges = test_utils::compute_subedges(ar_complex, repo, plane_index, orientations);
        auto r = ar_cut_2_face(ar_complex, 0, plane_index, orientations, subedges);
        REQUIRE(r[0] == 0);
        REQUIRE(r[1] == INVALID);
        REQUIRE(r[2] != INVALID);
    }
    SECTION("Tagent case: face")
    {
        // Inserting plane 5.
        size_t plane_index = 5;
        auto orientations = test_utils::compute_orientations(ar_complex, repo, plane_index);
        auto subedges = test_utils::compute_subedges(ar_complex, repo, plane_index, orientations);
        auto r = ar_cut_2_face(ar_complex, 0, plane_index, orientations, subedges);
        REQUIRE(r[0] == INVALID);
        REQUIRE(r[1] == INVALID);
        REQUIRE(r[2] == INVALID);
    }
    SECTION("No cut")
    {
        // Inserting plane 9.
        size_t plane_index = 9;
        auto orientations = test_utils::compute_orientations(ar_complex, repo, plane_index);
        auto subedges = test_utils::compute_subedges(ar_complex, repo, plane_index, orientations);
        auto r = ar_cut_2_face(ar_complex, 0, plane_index, orientations, subedges);
        REQUIRE(r[0] == 0);
        REQUIRE(r[1] == INVALID);
        REQUIRE(r[2] == INVALID);
    }
    SECTION("Cross cut")
    {
        // Inserting plane 8.
        size_t plane_index = 8;
        auto orientations = test_utils::compute_orientations(ar_complex, repo, plane_index);
        auto subedges = test_utils::compute_subedges(ar_complex, repo, plane_index, orientations);
        auto r = ar_cut_2_face(ar_complex, 0, plane_index, orientations, subedges);
        REQUIRE(r[0] != INVALID);
        REQUIRE(r[1] != INVALID);
        REQUIRE(r[2] != INVALID);

        const auto& f = ar_complex.faces[0];
        const auto& f0 = ar_complex.faces[r[0]];
        const auto& f1 = ar_complex.faces[r[1]];
        const auto& e = ar_complex.edges[r[2]];

        REQUIRE(f.supporting_plane == f0.supporting_plane);
        REQUIRE(f.supporting_plane == f1.supporting_plane);
        REQUIRE(f0.edges.size() == 3);
        REQUIRE(f1.edges.size() == 4);
        REQUIRE((e.supporting_planes[0] == plane_index || e.supporting_planes[1] == plane_index));

        check_face_orientation(ar_complex, r[0], r[1], r[2]);
    }
}
} // namespace

TEST_CASE("ar_cut_2_face", "[plane_interface]")
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
