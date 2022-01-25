#include <simplicial_arrangement/lookup_table.h>
#include <simplicial_arrangement/simplicial_arrangement.h>

#include <implicit_predicates/implicit_predicates.h>

#include <spdlog/spdlog.h>
#include <catch2/catch.hpp>

#include <cmath>

namespace {

template <typename Scalar, int DIM>
void validate_arrangement(simplicial_arrangement::Arrangement<DIM>& r,
    const std::vector<simplicial_arrangement::Plane<Scalar, DIM>>& planes)
{
    using namespace simplicial_arrangement;

    auto get_plane = [&](size_t i) {
        Plane<Scalar, DIM> P;
        if (i <= DIM) {
            for (size_t j = 0; j <= DIM; j++) {
                P[j] = (i == j) ? 1 : 0;
            }
        } else {
            P = planes[i - DIM - 1];
        }
        return P;
    };

    auto dot = [&](size_t pi, size_t pj) {
        auto Pi = get_plane(pi);
        auto Pj = get_plane(pj);
        Scalar s = 0;
        for (size_t i = 0; i <= DIM; i++) {
            s += Pi[i] * Pj[i];
        }
        return s;
    };

    size_t num_unique_planes = r.unique_planes.size();
    for (size_t i = 0; i < num_unique_planes; i++) {
        const auto& plane_group = r.unique_planes[i];
        if (plane_group.size() == 1) continue;
        const size_t size = plane_group.size();
        const auto& orientations = r.unique_plane_orientations;
        for (size_t j = 1; j < size; j++) {
            REQUIRE(
                r.unique_plane_indices[plane_group[0]] == r.unique_plane_indices[plane_group[j]]);
            if (orientations[plane_group[0]] == orientations[plane_group[j]]) {
                REQUIRE(dot(plane_group[0], plane_group[j]) > 0);
            } else {
                REQUIRE(dot(plane_group[0], plane_group[j]) < 0);
            }
        }
    }

    size_t num_faces = r.faces.size();
    for (size_t i = 0; i < num_faces; i++) {
        const auto& face = r.faces[i];
        if constexpr (DIM == 2) {
            REQUIRE(face.vertices.size() == 2);
        } else {
            REQUIRE(face.vertices.size() >= 3);
        }
        REQUIRE(face.supporting_plane != Arrangement<DIM>::None);

        if (!r.unique_planes.empty()) {
            REQUIRE(face.supporting_plane < r.unique_plane_indices.size());
            const size_t unique_plane_index = r.unique_plane_indices[face.supporting_plane];
            REQUIRE(unique_plane_index < r.unique_planes.size());
            const auto& supporting_planes = r.unique_planes[unique_plane_index];
            REQUIRE(supporting_planes.size() >= 1);

            auto itr = std::find(
                supporting_planes.begin(), supporting_planes.end(), face.supporting_plane);
            REQUIRE(itr != supporting_planes.end());
        }

        for (size_t vid : face.vertices) {
            REQUIRE(vid < r.vertices.size());
        }
    }

    size_t num_cells = r.cells.size();
    for (size_t i = 0; i < num_cells; i++) {
        const auto& cell = r.cells[i];
        REQUIRE(cell.faces.size() >= DIM + 1);
        size_t num_faces = cell.faces.size();
        for (size_t j = 0; j < num_faces; j++) {
            const auto& face = r.faces[cell.faces[j]];

            if (face.positive_cell == Arrangement<DIM>::None ||
                face.negative_cell == Arrangement<DIM>::None) {
                size_t plane_id = DIM + 1;
                if (!r.unique_planes.empty()) {
                    const auto& supporting_planes =
                        r.unique_planes[r.unique_plane_indices[face.supporting_plane]];
                    auto itr = std::min_element(supporting_planes.begin(), supporting_planes.end());
                    plane_id = *itr;
                } else {
                    plane_id = face.supporting_plane;
                }
                REQUIRE(plane_id <= DIM);

                // Cell must be on the positive side of the tet boundary planes.
                if (plane_id != face.supporting_plane) {
                    if (r.unique_plane_orientations[face.supporting_plane] ==
                        r.unique_plane_orientations[plane_id]) {
                        REQUIRE(face.positive_cell == i);
                    } else {
                        REQUIRE(face.negative_cell == i);
                    }
                }
            }
        }
    }
}

template <int DIM>
void assert_equivalent(const simplicial_arrangement::Arrangement<DIM>& r1,
    const simplicial_arrangement::Arrangement<DIM>& r2)
{
    using namespace simplicial_arrangement;

    REQUIRE(r1.vertices.size() == r2.vertices.size());
    REQUIRE(r1.faces.size() == r2.faces.size());
    REQUIRE(r1.cells.size() == r2.cells.size());

    bool r1_all_distinct = r1.unique_planes.size() == r1.unique_plane_indices.size();
    bool r2_all_distinct = r2.unique_planes.size() == r2.unique_plane_indices.size();
    REQUIRE(r1_all_distinct == r2_all_distinct);

    if (!r1_all_distinct) {
        REQUIRE(r1.unique_planes.size() == r2.unique_planes.size());

        const size_t num_unique_planes = r1.unique_planes.size();
        std::vector<size_t> r1_support_face_count(num_unique_planes, 0);
        std::vector<size_t> r2_support_face_count(num_unique_planes, 0);

        for (const auto& f : r1.faces) {
            size_t uid = r1.unique_plane_indices[f.supporting_plane];
            r1_support_face_count[uid]++;
        }
        for (const auto& f : r2.faces) {
            size_t uid = r2.unique_plane_indices[f.supporting_plane];
            r2_support_face_count[uid]++;
        }

        const size_t num_planes = r1.unique_plane_indices.size();
        for (size_t i = 0; i < num_planes; i++) {
            size_t r1_uid = r1.unique_plane_indices[i];
            size_t r2_uid = r2.unique_plane_indices[i];
            CAPTURE(r1.unique_plane_indices,
                r2.unique_plane_indices,
                r1_support_face_count,
                r2_support_face_count);
            CHECK(r1.unique_planes[r1_uid].size() == r1.unique_planes[r2_uid].size());
            CHECK(r1_support_face_count[r1_uid] == r2_support_face_count[r2_uid]);
        }
    }
}

template <typename Scalar>
void test_2D()
{
    using namespace simplicial_arrangement;

    auto count_num_halfedges = [](Arrangement<2>& arrangement) {
        size_t count = 0;
        for (const auto& cell : arrangement.cells) {
            count += cell.faces.size();
        }
        return count;
    };

    std::vector<Plane<Scalar, 2>> planes;

    SECTION("0 implicit")
    {
        auto arrangement = compute_arrangement(planes);
        REQUIRE(arrangement.cells.size() == 1);
        REQUIRE(count_num_halfedges(arrangement) == 3);
        REQUIRE(arrangement.vertices.size() == 3);
        validate_arrangement(arrangement, planes);
    }

    SECTION("1 implicit")
    {
        SECTION("normal case")
        {
            planes.push_back({1, -1, -1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(arrangement.cells.size() == 2);
            REQUIRE(count_num_halfedges(arrangement) == 7);
            REQUIRE(arrangement.vertices.size() == 5);
            validate_arrangement(arrangement, planes);
        }
        SECTION("passing through one corner")
        {
            planes.push_back({1, -1, 0});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(arrangement.cells.size() == 2);
            REQUIRE(count_num_halfedges(arrangement) == 6);
            REQUIRE(arrangement.vertices.size() == 4);
            validate_arrangement(arrangement, planes);
        }
        SECTION("passing through two corners")
        {
            planes.push_back({1, 0, 0});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(arrangement.cells.size() == 1);
            REQUIRE(count_num_halfedges(arrangement) == 3);
            REQUIRE(arrangement.vertices.size() == 3);
            REQUIRE(arrangement.unique_planes.size() == 3);
            validate_arrangement(arrangement, planes);
        }
    }

    SECTION("2 implicits")
    {
        SECTION("not intersecting")
        {
            planes.push_back({1, -1, -1});
            planes.push_back({2, -1, -1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(arrangement.cells.size() == 3);
            REQUIRE(count_num_halfedges(arrangement) == 11);
            REQUIRE(arrangement.vertices.size() == 7);
            validate_arrangement(arrangement, planes);
        }

        SECTION("crossing")
        {
            planes.push_back({1, -1, -1});
            planes.push_back({-1, 2, -1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(arrangement.cells.size() == 4);
            REQUIRE(count_num_halfedges(arrangement) == 15);
            REQUIRE(arrangement.vertices.size() == 8);
            validate_arrangement(arrangement, planes);
        }

        SECTION("touching")
        {
            planes.push_back({1, -1, -1});
            planes.push_back({-1, 1, -1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(arrangement.cells.size() == 3);
            REQUIRE(count_num_halfedges(arrangement) == 10);
            REQUIRE(arrangement.vertices.size() == 6);
            validate_arrangement(arrangement, planes);
        }

        SECTION("overlapping")
        {
            planes.push_back({1, -1, -1});
            planes.push_back({1, -1, -1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(arrangement.cells.size() == 2);
            REQUIRE(count_num_halfedges(arrangement) == 7);
            REQUIRE(arrangement.vertices.size() == 5);
            REQUIRE(arrangement.unique_planes.size() == 4);
            validate_arrangement(arrangement, planes);
        }
    }
    SECTION("3 implicits")
    {
        SECTION("non-intersecting")
        {
            planes.push_back({1, -1, -1});
            planes.push_back({2, -1, -1});
            planes.push_back({3, -1, -1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(arrangement.cells.size() == 4);
            REQUIRE(count_num_halfedges(arrangement) == 15);
            REQUIRE(arrangement.vertices.size() == 9);
            validate_arrangement(arrangement, planes);
        }

        SECTION("intersecting")
        {
            SECTION("case 0")
            {
                planes.push_back({3, -1, -1});
                planes.push_back({-1, 3, -1});
                planes.push_back({-1, -1, 3});
                auto arrangement = compute_arrangement(planes);
                REQUIRE(arrangement.cells.size() == 7);
                REQUIRE(count_num_halfedges(arrangement) == 27);
                REQUIRE(arrangement.vertices.size() == 12);
                validate_arrangement(arrangement, planes);
            }

            SECTION("case 1")
            {
                planes.push_back({2, -1, -1});
                planes.push_back({-1, 2, -1});
                planes.push_back({-1, -1, 2});
                auto arrangement = compute_arrangement(planes);
                REQUIRE(arrangement.cells.size() == 6);
                REQUIRE(count_num_halfedges(arrangement) == 21);
                REQUIRE(arrangement.vertices.size() == 10);
                validate_arrangement(arrangement, planes);
            }

            SECTION("case 2")
            {
                planes.push_back({1, -1, -1});
                planes.push_back({-1, 1, -1});
                planes.push_back({-1, -1, 1});
                auto arrangement = compute_arrangement(planes);
                REQUIRE(arrangement.cells.size() == 4);
                REQUIRE(count_num_halfedges(arrangement) == 12);
                REQUIRE(arrangement.vertices.size() == 6);
                validate_arrangement(arrangement, planes);
            }

            SECTION("case 3")
            {
                planes.push_back({1, -2, -2});
                planes.push_back({-2, 1, -2});
                planes.push_back({-2, -2, 1});
                auto arrangement = compute_arrangement(planes);
                REQUIRE(arrangement.cells.size() == 4);
                REQUIRE(count_num_halfedges(arrangement) == 15);
                REQUIRE(arrangement.vertices.size() == 9);
                validate_arrangement(arrangement, planes);
            }

            SECTION("case 4")
            {
                planes.push_back({3, -1, -1});
                planes.push_back({-1, 3, -1});
                planes.push_back({-1, -1, 2});
                auto arrangement = compute_arrangement(planes);
                REQUIRE(arrangement.cells.size() == 7);
                REQUIRE(count_num_halfedges(arrangement) == 27);
                REQUIRE(arrangement.vertices.size() == 12);
                validate_arrangement(arrangement, planes);
            }

            SECTION("case 5")
            {
                planes.push_back({3, -1, -1});
                planes.push_back({-3, 1, 1});
                planes.push_back({-1, -1, 2});
                auto arrangement = compute_arrangement(planes);
                REQUIRE(arrangement.cells.size() == 4);
                REQUIRE(count_num_halfedges(arrangement) == 15);
                REQUIRE(arrangement.vertices.size() == 8);
                REQUIRE(arrangement.unique_planes.size() == 5);
                validate_arrangement(arrangement, planes);
            }
        }
    }
}

template <typename Scalar>
void test_3D()
{
    using namespace simplicial_arrangement;

    auto count_num_cells = [](auto& arrangement) { return arrangement.cells.size(); };

    auto count_num_half_faces = [](auto& arrangement) {
        size_t count = 0;
        for (const auto& cell : arrangement.cells) {
            count += cell.faces.size();
        }
        return count;
    };

    auto count_num_half_edges = [](auto& arrangement) {
        size_t count = 0;
        for (const auto& cell : arrangement.cells) {
            for (const auto fid : cell.faces) {
                count += arrangement.faces[fid].vertices.size();
            }
        }
        return count;
    };

    std::vector<Plane<Scalar, 3>> planes;

    SECTION("0 implicit")
    {
        auto arrangement = compute_arrangement(planes);
        REQUIRE(count_num_cells(arrangement) == 1);
        REQUIRE(count_num_half_faces(arrangement) == 4);
        REQUIRE(count_num_half_edges(arrangement) == 12);
        REQUIRE(arrangement.vertices.size() == 4);
        validate_arrangement(arrangement, planes);
    }

    SECTION("1 implicit")
    {
        SECTION("quad cross section")
        {
            planes.push_back({1, 1, -1, -1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 2);
            REQUIRE(count_num_half_faces(arrangement) == 10);
            REQUIRE(count_num_half_edges(arrangement) == 36);
            REQUIRE(arrangement.vertices.size() == 8);
            validate_arrangement(arrangement, planes);
        }
        SECTION("tri cross section")
        {
            planes.push_back({1, -1, -1, -1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 2);
            REQUIRE(count_num_half_faces(arrangement) == 9);
            REQUIRE(count_num_half_edges(arrangement) == 30);
            REQUIRE(arrangement.vertices.size() == 7);
            validate_arrangement(arrangement, planes);
        }
        SECTION("cut through a vertex")
        {
            planes.push_back({0, -1, -1, 1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 2);
            REQUIRE(count_num_half_faces(arrangement) == 9);
            REQUIRE(count_num_half_edges(arrangement) == 28);
            REQUIRE(arrangement.vertices.size() == 6);
            validate_arrangement(arrangement, planes);
        }
        SECTION("cut through an edge")
        {
            planes.push_back({0, 0, -1, 1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 2);
            REQUIRE(count_num_half_faces(arrangement) == 8);
            REQUIRE(count_num_half_edges(arrangement) == 24);
            REQUIRE(arrangement.vertices.size() == 5);
            validate_arrangement(arrangement, planes);
        }
        SECTION("cut through a face")
        {
            planes.push_back({0, 0, 0, 1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 1);
            REQUIRE(count_num_half_faces(arrangement) == 4);
            REQUIRE(count_num_half_edges(arrangement) == 12);
            REQUIRE(arrangement.vertices.size() == 4);
            REQUIRE(arrangement.unique_planes.size() == 4);
            validate_arrangement(arrangement, planes);
        }
    }

    SECTION("2 implicits")
    {
        SECTION("not intersecting cuts")
        {
            planes.push_back({1, -2, -2, -2});
            planes.push_back({-2, 1, -2, -2});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 3);
            REQUIRE(count_num_half_faces(arrangement) == 14);
            REQUIRE(count_num_half_edges(arrangement) == 48);
            REQUIRE(arrangement.vertices.size() == 10);
            validate_arrangement(arrangement, planes);
        }
        SECTION("crossing")
        {
            planes.push_back({2, -1, -1, -1});
            planes.push_back({-1, 2, -1, -1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 4);
            REQUIRE(count_num_half_faces(arrangement) == 20);
            REQUIRE(count_num_half_edges(arrangement) == 72);
            REQUIRE(arrangement.vertices.size() == 12);
            validate_arrangement(arrangement, planes);
        }
        SECTION("touching")
        {
            planes.push_back({1, -1, -1, -1});
            planes.push_back({-1, 1, -1, -1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 3);
            REQUIRE(count_num_half_faces(arrangement) == 14);
            REQUIRE(count_num_half_edges(arrangement) == 46);
            REQUIRE(arrangement.vertices.size() == 9);
            validate_arrangement(arrangement, planes);
        }
        SECTION("duplicated cuts")
        {
            planes.push_back({1, 1, -1, -1});
            planes.push_back({2, 2, -2, -2});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 2);
            REQUIRE(count_num_half_faces(arrangement) == 10);
            REQUIRE(count_num_half_edges(arrangement) == 36);
            REQUIRE(arrangement.vertices.size() == 8);
            REQUIRE(arrangement.unique_planes.size() == 5);
            validate_arrangement(arrangement, planes);
        }
        SECTION("2 cuts through the same vertex")
        {
            planes.push_back({0, 0, 1, -1});
            planes.push_back({0, 2, -2, 0});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 4);
            REQUIRE(count_num_half_faces(arrangement) == 17);
            REQUIRE(count_num_half_edges(arrangement) == 52);
            REQUIRE(arrangement.vertices.size() == 7);
            validate_arrangement(arrangement, planes);
        }
        SECTION("2 cuts through faces")
        {
            planes.push_back({0, 0, 0, 1});
            planes.push_back({0, 0, 0, 2});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 1);
            REQUIRE(count_num_half_faces(arrangement) == 4);
            REQUIRE(count_num_half_edges(arrangement) == 12);
            REQUIRE(arrangement.vertices.size() == 4);
            REQUIRE(arrangement.unique_planes.size() == 4);
            validate_arrangement(arrangement, planes);
        }
    }
    SECTION("3 implicits")
    {
        SECTION("3 duplicated cuts")
        {
            planes.push_back({1, 1, 1, -1});
            planes.push_back({3, 3, 3, -3});
            planes.push_back({-3, -3, -3, 3});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 2);
            REQUIRE(count_num_half_faces(arrangement) == 9);
            REQUIRE(count_num_half_edges(arrangement) == 30);
            REQUIRE(arrangement.vertices.size() == 7);
            REQUIRE(arrangement.unique_planes.size() == 5);
            validate_arrangement(arrangement, planes);
        }
        SECTION("3 tangent cuts")
        {
            planes.push_back({1, 0, 0, 0});
            planes.push_back({0, 0, 0, -3});
            planes.push_back({0, 0, -30, 0});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 1);
            REQUIRE(count_num_half_faces(arrangement) == 4);
            REQUIRE(count_num_half_edges(arrangement) == 12);
            REQUIRE(arrangement.vertices.size() == 4);
            REQUIRE(arrangement.unique_planes.size() == 4);
            validate_arrangement(arrangement, planes);
        }
        SECTION("3 parallel cuts")
        {
            planes.push_back({-1, 4, 4, 4});
            planes.push_back({-2, 4, 4, 4});
            planes.push_back({-3, 4, 4, 4});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 4);
            REQUIRE(count_num_half_faces(arrangement) == 19);
            REQUIRE(count_num_half_edges(arrangement) == 66);
            REQUIRE(arrangement.vertices.size() == 13);
            validate_arrangement(arrangement, planes);
        }
        SECTION("3 cuts intersect at a line")
        {
            planes.push_back({-1, 1, 0, 0});
            planes.push_back({-2, 0, 2, 0});
            planes.push_back({0, -4, 4, 0});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 6);
            REQUIRE(count_num_half_faces(arrangement) == 24);
            REQUIRE(count_num_half_edges(arrangement) == 72);
            REQUIRE(arrangement.vertices.size() == 8);
            validate_arrangement(arrangement, planes);
        }
        SECTION("3 cuts non-intersecting")
        {
            planes.push_back({1, -10, -10, -10});
            planes.push_back({-10, 1, -10, -10});
            planes.push_back({-10, -10, 1, -10});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 4);
            REQUIRE(count_num_half_faces(arrangement) == 19);
            REQUIRE(count_num_half_edges(arrangement) == 66);
            REQUIRE(arrangement.vertices.size() == 13);
            validate_arrangement(arrangement, planes);
        }
        SECTION("3 cuts intersect at a point")
        {
            planes.push_back({10, -1, -1, -1});
            planes.push_back({-1, 10, -1, -1});
            planes.push_back({-1, -1, 10, -1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 8);
            REQUIRE(count_num_half_faces(arrangement) == 43);
            REQUIRE(count_num_half_edges(arrangement) == 162);
            REQUIRE(arrangement.vertices.size() == 20);
            validate_arrangement(arrangement, planes);
        }
        SECTION("3 cuts intersect at a vertex")
        {
            planes.push_back({1, -1, -1, 0});
            planes.push_back({-1, 1, -1, 0});
            planes.push_back({-1, -1, 1, 0});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 4);
            REQUIRE(count_num_half_faces(arrangement) == 16);
            REQUIRE(count_num_half_edges(arrangement) == 48);
            REQUIRE(arrangement.vertices.size() == 7);
            validate_arrangement(arrangement, planes);
        }
        SECTION("3 cuts intersect at an edge")
        {
            planes.push_back({1, -1, -1, -1});
            planes.push_back({-1, 1, -2, -2});
            planes.push_back({1, -1, 3, 3});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 4);
            REQUIRE(count_num_half_faces(arrangement) == 19);
            REQUIRE(count_num_half_edges(arrangement) == 62);
            REQUIRE(arrangement.vertices.size() == 11);
            validate_arrangement(arrangement, planes);
        }
    }
}

} // namespace

TEST_CASE("Arrangement 2D", "[arrangement][2D]")
{
    using namespace simplicial_arrangement;
    disable_lookup_table();
#ifndef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
    SECTION("Int") { test_2D<Int>(); }
#endif
    SECTION("double") { test_2D<double>(); }
}

TEST_CASE("Arrangement 3D", "[arrangement][3d]")
{
    using namespace simplicial_arrangement;
    disable_lookup_table();
#ifndef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
    SECTION("Int") { test_3D<Int>(); }
#endif
    SECTION("double") { test_3D<double>(); }
}

TEST_CASE("Lookup 3D", "[lookup][3D]")
{
    using namespace simplicial_arrangement;
    bool loaded = load_lookup_table(ARRANGEMENT);
    REQUIRE(loaded == true);

#ifndef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
    SECTION("Int") { test_3D<Int>(); }
#endif
    SECTION("double") { test_3D<double>(); }

    SECTION("Consistency check")
    {
        SECTION("1 cut")
        {
            using Scalar = double;
            constexpr Scalar EPS = std::numeric_limits<Scalar>::epsilon();
            std::vector<Plane<Scalar, 3>> planes;

            SECTION("Case 1")
            {
                planes.push_back({EPS, -EPS, EPS, -EPS});

                enable_lookup_table();
                auto r1 = compute_arrangement(planes);
                disable_lookup_table();
                auto r2 = compute_arrangement(planes);

                assert_equivalent(r1, r2);
            }
            SECTION("Case 2")
            {
                planes.push_back({EPS, EPS, EPS, -EPS});

                enable_lookup_table();
                auto r1 = compute_arrangement(planes);
                disable_lookup_table();
                auto r2 = compute_arrangement(planes);

                assert_equivalent(r1, r2);
            }
            SECTION("Case 3")
            {
                planes.push_back({EPS, EPS, EPS, 0});

                enable_lookup_table();
                auto r1 = compute_arrangement(planes);
                disable_lookup_table();
                auto r2 = compute_arrangement(planes);

                assert_equivalent(r1, r2);
            }
            SECTION("Case 4")
            {
                planes.push_back({EPS, -EPS, 0, 0});

                enable_lookup_table();
                auto r1 = compute_arrangement(planes);
                disable_lookup_table();
                auto r2 = compute_arrangement(planes);

                assert_equivalent(r1, r2);
            }
            SECTION("Case 5")
            {
                planes.push_back({EPS, 0, 0, 0});

                enable_lookup_table();
                auto r1 = compute_arrangement(planes);
                disable_lookup_table();
                auto r2 = compute_arrangement(planes);

                assert_equivalent(r1, r2);
            }
            SECTION("Case 6")
            {
                planes.push_back({EPS, EPS, EPS, EPS});

                enable_lookup_table();
                auto r1 = compute_arrangement(planes);
                disable_lookup_table();
                auto r2 = compute_arrangement(planes);

                assert_equivalent(r1, r2);
            }
        }
        SECTION("2 cut")
        {
            using Scalar = double;
            constexpr Scalar EPS = std::numeric_limits<Scalar>::epsilon();
            constexpr Scalar l = 1e12;
            std::vector<Plane<Scalar, 3>> planes;

            SECTION("Case 1")
            {
                planes.push_back({EPS, EPS, EPS, l});
                planes.push_back({-EPS, -EPS, -EPS, -l});

                enable_lookup_table();
                auto r1 = compute_arrangement(planes);
                disable_lookup_table();
                auto r2 = compute_arrangement(planes);

                assert_equivalent(r1, r2);
            }
            SECTION("Case 2")
            {
                planes.push_back({EPS, -EPS, 0, 0});
                planes.push_back({-l, l, 0, 0});

                enable_lookup_table();
                auto r1 = compute_arrangement(planes);
                disable_lookup_table();
                auto r2 = compute_arrangement(planes);

                assert_equivalent(r1, r2);
            }
            SECTION("Case 3")
            {
                planes.push_back({EPS, -EPS, 1, 0});
                planes.push_back({0, l, 2 * l, -l});

                enable_lookup_table();
                auto r1 = compute_arrangement(planes);
                disable_lookup_table();
                auto r2 = compute_arrangement(planes);

                assert_equivalent(r1, r2);
            }
            SECTION("Case 4")
            {
                planes.push_back({EPS, EPS, EPS, -l});
                planes.push_back({-l, -l, -l, EPS});

                enable_lookup_table();
                auto r1 = compute_arrangement(planes);
                disable_lookup_table();
                auto r2 = compute_arrangement(planes);

                assert_equivalent(r1, r2);
            }
        }
    }
}

/**
 * Note that the following tests will fail when
 * `SIMPLICIAL_ARRANGMEENT_NON_ROBUST` compiler flag is set.
 */
TEST_CASE("Arrangement degeneracy", "[ar][3D]")
{
    using namespace simplicial_arrangement;
    disable_lookup_table();

    using Scalar = double;
    constexpr Scalar EPS = std::numeric_limits<Scalar>::epsilon();
    constexpr Scalar LARGE = M_PI * 1e6;
    std::vector<Plane<Scalar, 3>> planes;

    SECTION("Near parallel planes")
    {
        // Reference case:
        planes.push_back({-1.1, -1.2, 1.3, 1.4});
        planes.push_back({-1.2, -1.1, 1.3, 1.4});
        planes.push_back({-1.3, -1.2, 1.1, 1.4});
        planes.push_back({-1.4, -1.2, 1.3, 1.1});
        const auto r1 = compute_arrangement(planes);

        SECTION("With normal values")
        {
            Scalar v1 = M_PI;
            Scalar v2 = M_PI + 1e-4;
            Scalar v3 = M_PI + 1e-3;
            Scalar v4 = M_PI + 2e-3;
            planes.clear();
            planes.push_back({-v1, -v2, v3, v4});
            planes.push_back({-v2, -v1, v3, v4});
            planes.push_back({-v3, -v2, v1, v4});
            planes.push_back({-v4, -v2, v3, v1});
            const auto r2 = compute_arrangement(planes);
            assert_equivalent(r1, r2);
        }

        SECTION("With small values")
        {
            const Scalar EPS2 = std::nextafter(EPS, 1.0);
            const Scalar EPS3 = std::nextafter(EPS2, 1.0);
            const Scalar EPS4 = std::nextafter(EPS3, 1.0);
            planes.clear();
            planes.push_back({-EPS, -EPS2, EPS3, EPS4});
            planes.push_back({-EPS2, -EPS, EPS3, EPS4});
            planes.push_back({-EPS3, -EPS2, EPS, EPS4});
            planes.push_back({-EPS4, -EPS2, EPS3, EPS});
            // spdlog::set_level(spdlog::level::debug);
            const auto r2 = compute_arrangement(planes);
            assert_equivalent(r1, r2);
        }

        SECTION("With large values")
        {
            const Scalar LARGE2 = std::nextafter(LARGE, 1.0);
            const Scalar LARGE3 = std::nextafter(LARGE2, 1.0);
            const Scalar LARGE4 = std::nextafter(LARGE3, 1.0);
            planes.clear();
            planes.push_back({-LARGE, -LARGE2, LARGE3, LARGE4});
            planes.push_back({-LARGE2, -LARGE, LARGE3, LARGE4});
            planes.push_back({-LARGE3, -LARGE2, LARGE, LARGE4});
            planes.push_back({-LARGE4, -LARGE2, LARGE3, LARGE});
            const auto r3 = compute_arrangement(planes);
            assert_equivalent(r1, r3);
        }
    }

    SECTION("4 planes intersecting at the same point")
    {
        // Reference case:
        planes.push_back({-1, -1, 1, 1});
        planes.push_back({-1, 1, -1, 1});
        planes.push_back({-1, 1, 1, -1});
        planes.push_back({2, -1, -2, 1});
        const auto r1 = compute_arrangement(planes);

        SECTION("With normal values")
        {
            const Scalar v1 = M_PI;
            const Scalar v2 = M_PI / 3;
            const Scalar v3 = 1.1;
            const Scalar v4 = std::sqrt(99);
            planes.clear();
            planes.push_back({-v1, -v1, v1, v1});
            planes.push_back({-v2, v2, -v2, v2});
            planes.push_back({-v3, v3, v3, -v3});
            planes.push_back({v4, -v2, -v4, v2});
            const auto r2 = compute_arrangement(planes);
            assert_equivalent(r1, r2);
        }
    }

    SECTION("4 planes in different order")
    {
        SECTION("Point intersection")
        {
            planes.push_back(
                {6.926218366896812, -3.734529738225069, 0.4909632573147267, -3.68265188598647});
            planes.push_back(
                {-2.241785194788779, 3.394920808940942, 8.710781454093603, -9.863917068245765});
            planes.push_back(
                {9.980810309305472, -5.278220474036615, -2.068385476747814, -2.634204358521043});
            planes.push_back(
                {9.94369616460531, 8.65114722736331, -7.437511044553879, -11.15733234741474});

            const auto r1 = compute_arrangement(planes);
            std::reverse(planes.begin(), planes.end());
            const auto r2 = compute_arrangement(planes);
            assert_equivalent(r1, r2);
        }
        SECTION("Segment intersection")
        {
            planes.push_back(
                {17.24675904665293, -26.21427764215589, 0.6882781443530117, 8.279240451149956});
            planes.push_back(
                {-11.94785331394204, 18.48725103287373, -1.130942123921345, -5.408455595010345});
            planes.push_back(
                {-2.752603223595615, 5.996169704505958, -3.734529738225069, 0.4909632573147267});
            planes.push_back(
                {22.56321818788723, -38.20021800887764, 8.710781454093603, 6.926218366896812});

            const auto r1 = compute_arrangement(planes);
            std::reverse(planes.begin(), planes.end());
            const auto r2 = compute_arrangement(planes);
            assert_equivalent(r1, r2);
        }
        SECTION("Tri intersection")
        {
            planes.push_back(
                {-6.76762939664222, 6.767629402413544, 6.767629396432295, 6.767629409302009});
            planes.push_back(
                {-4.112810564622128, 4.11281055783204, 4.112810554094155, 4.11281055555015});
            planes.push_back(
                {5.639169256076264, -5.639169262552152, -5.639169251034964, -5.639169260805547});
            planes.push_back(
                {-1.432351809908365, 1.432351800879061, 1.432351815479439, 1.432351815092373});

            const auto r1 = compute_arrangement(planes);
            std::reverse(planes.begin(), planes.end());
            const auto r2 = compute_arrangement(planes);
            assert_equivalent(r1, r2);
        }
        SECTION("Quad intersection")
        {
            planes.push_back(
                {9.406140896875673, 9.406140889947853, -9.406140902713707, -9.40614090529332});
            planes.push_back(
                {-9.912418053618833, -9.912418066025847, 9.912418049659827, 9.912418049968453});
            planes.push_back(
                {0.3326368104448145, 0.3326368139110936, -0.3326368178811913, -0.3326368006459226});
            planes.push_back(
                {-7.077317534324832, -7.077317531984346, 7.077317548576322, 7.07731754191976});

            const auto r1 = compute_arrangement(planes);
            std::reverse(planes.begin(), planes.end());
            const auto r2 = compute_arrangement(planes);
            assert_equivalent(r1, r2);
        }
    }
}
