#include <catch2/catch.hpp>

#include <simplicial_arrangement/simplicial_arrangement.h>

namespace {

template <int DIM>
void validate_arrangement(simplicial_arrangement::Arrangement<DIM>& r)
{
    using namespace simplicial_arrangement;

    size_t num_faces = r.faces.size();
    for (size_t i = 0; i < num_faces; i++) {
        const auto& face = r.faces[i];
        if constexpr (DIM == 2) {
            REQUIRE(face.vertices.size() == 2);
        } else {
            REQUIRE(face.vertices.size() >= 3);
        }
        REQUIRE(face.supporting_plane != Arrangement<DIM>::None);
        REQUIRE(face.supporting_plane < r.unique_plane_indices.size());

        const size_t unique_plane_index = r.unique_plane_indices[face.supporting_plane];
        REQUIRE(unique_plane_index < r.unique_planes.size());
        const auto& supporting_planes = r.unique_planes[unique_plane_index];
        REQUIRE(supporting_planes.size() >= 1);

        auto itr =
            std::find(supporting_planes.begin(), supporting_planes.end(), face.supporting_plane);
        REQUIRE(itr != supporting_planes.end());

        for (size_t vid : face.vertices) {
            REQUIRE(vid < r.vertices.size());
        }
    }

    size_t num_cells = r.cells.size();
    for (size_t i = 0; i < num_cells; i++) {
        const auto& cell = r.cells[i];
        REQUIRE(cell.faces.size() >= DIM + 1);
        REQUIRE(cell.faces.size() == cell.face_orientations.size());
        size_t num_faces = cell.faces.size();
        for (size_t j = 0; j < num_faces; j++) {
            const auto& face = r.faces[cell.faces[j]];
            if (cell.face_orientations[j]) {
                REQUIRE(face.positive_cell == i);
            } else {
                REQUIRE(face.negative_cell == i);
            }

            if (face.positive_cell == Arrangement<DIM>::None ||
                face.negative_cell == Arrangement<DIM>::None) {
                const auto& supporting_planes =
                    r.unique_planes[r.unique_plane_indices[face.supporting_plane]];
                auto itr = std::min_element(supporting_planes.begin(), supporting_planes.end());
                REQUIRE(*itr <= DIM);

                size_t k = itr - supporting_planes.begin();
                itr = std::find(
                    supporting_planes.begin(), supporting_planes.end(), face.supporting_plane);
                REQUIRE(itr != supporting_planes.end());
                size_t l = itr - supporting_planes.begin();

                if (r.unique_plane_orientations[k] == r.unique_plane_orientations[l]) {
                    REQUIRE(cell.face_orientations[j]);
                } else {
                    REQUIRE(!cell.face_orientations[j]);
                }
            }
        }
    }
};

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
        REQUIRE(arrangement.unique_planes.size() == 3);
        validate_arrangement(arrangement);
    }

    SECTION("1 implicit")
    {
        SECTION("normal case")
        {
            planes.push_back({1, -1, -1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(arrangement.unique_plane_indices.size() == 4);
            REQUIRE(arrangement.cells.size() == 2);
            REQUIRE(count_num_halfedges(arrangement) == 7);
            REQUIRE(arrangement.vertices.size() == 5);
            REQUIRE(arrangement.unique_planes.size() == 4);
            validate_arrangement(arrangement);
        }
        SECTION("passing through one corner")
        {
            planes.push_back({1, -1, 0});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(arrangement.unique_plane_indices.size() == 4);
            REQUIRE(arrangement.cells.size() == 2);
            REQUIRE(count_num_halfedges(arrangement) == 6);
            REQUIRE(arrangement.vertices.size() == 4);
            REQUIRE(arrangement.unique_planes.size() == 4);
            validate_arrangement(arrangement);
        }
        SECTION("passing through two corners")
        {
            planes.push_back({1, 0, 0});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(arrangement.unique_plane_indices.size() == 4);
            REQUIRE(arrangement.cells.size() == 1);
            REQUIRE(count_num_halfedges(arrangement) == 3);
            REQUIRE(arrangement.vertices.size() == 3);
            REQUIRE(arrangement.unique_planes.size() == 3);
            validate_arrangement(arrangement);
        }
    }

    SECTION("2 implicits")
    {
        SECTION("not intersecting")
        {
            planes.push_back({1, -1, -1});
            planes.push_back({2, -1, -1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(arrangement.unique_plane_indices.size() == 5);
            REQUIRE(arrangement.cells.size() == 3);
            REQUIRE(count_num_halfedges(arrangement) == 11);
            REQUIRE(arrangement.vertices.size() == 7);
            REQUIRE(arrangement.unique_planes.size() == 5);
            validate_arrangement(arrangement);
        }

        SECTION("crossing")
        {
            planes.push_back({1, -1, -1});
            planes.push_back({-1, 2, -1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(arrangement.unique_plane_indices.size() == 5);
            REQUIRE(arrangement.cells.size() == 4);
            REQUIRE(count_num_halfedges(arrangement) == 15);
            REQUIRE(arrangement.vertices.size() == 8);
            REQUIRE(arrangement.unique_planes.size() == 5);
            validate_arrangement(arrangement);
        }

        SECTION("touching")
        {
            planes.push_back({1, -1, -1});
            planes.push_back({-1, 1, -1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(arrangement.unique_plane_indices.size() == 5);
            REQUIRE(arrangement.cells.size() == 3);
            REQUIRE(count_num_halfedges(arrangement) == 10);
            REQUIRE(arrangement.vertices.size() == 6);
            REQUIRE(arrangement.unique_planes.size() == 5);
            validate_arrangement(arrangement);
        }

        SECTION("overlapping")
        {
            planes.push_back({1, -1, -1});
            planes.push_back({1, -1, -1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(arrangement.unique_plane_indices.size() == 5);
            REQUIRE(arrangement.cells.size() == 2);
            REQUIRE(count_num_halfedges(arrangement) == 7);
            REQUIRE(arrangement.vertices.size() == 5);
            REQUIRE(arrangement.unique_planes.size() == 4);
            validate_arrangement(arrangement);
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
            REQUIRE(arrangement.unique_plane_indices.size() == 6);
            REQUIRE(arrangement.cells.size() == 4);
            REQUIRE(count_num_halfedges(arrangement) == 15);
            REQUIRE(arrangement.vertices.size() == 9);
            REQUIRE(arrangement.unique_planes.size() == 6);
            validate_arrangement(arrangement);
        }

        SECTION("intersecting")
        {
            SECTION("case 0")
            {
                planes.push_back({3, -1, -1});
                planes.push_back({-1, 3, -1});
                planes.push_back({-1, -1, 3});
                auto arrangement = compute_arrangement(planes);
                REQUIRE(arrangement.unique_plane_indices.size() == 6);
                REQUIRE(arrangement.cells.size() == 7);
                REQUIRE(count_num_halfedges(arrangement) == 27);
                REQUIRE(arrangement.vertices.size() == 12);
                REQUIRE(arrangement.unique_planes.size() == 6);
                validate_arrangement(arrangement);
            }

            SECTION("case 1")
            {
                planes.push_back({2, -1, -1});
                planes.push_back({-1, 2, -1});
                planes.push_back({-1, -1, 2});
                auto arrangement = compute_arrangement(planes);
                REQUIRE(arrangement.unique_plane_indices.size() == 6);
                REQUIRE(arrangement.cells.size() == 6);
                REQUIRE(count_num_halfedges(arrangement) == 21);
                REQUIRE(arrangement.vertices.size() == 10);
                REQUIRE(arrangement.unique_planes.size() == 6);
                validate_arrangement(arrangement);
            }

            SECTION("case 2")
            {
                planes.push_back({1, -1, -1});
                planes.push_back({-1, 1, -1});
                planes.push_back({-1, -1, 1});
                auto arrangement = compute_arrangement(planes);
                REQUIRE(arrangement.unique_plane_indices.size() == 6);
                REQUIRE(arrangement.cells.size() == 4);
                REQUIRE(count_num_halfedges(arrangement) == 12);
                REQUIRE(arrangement.vertices.size() == 6);
                REQUIRE(arrangement.unique_planes.size() == 6);
                validate_arrangement(arrangement);
            }

            SECTION("case 3")
            {
                planes.push_back({1, -2, -2});
                planes.push_back({-2, 1, -2});
                planes.push_back({-2, -2, 1});
                auto arrangement = compute_arrangement(planes);
                REQUIRE(arrangement.unique_plane_indices.size() == 6);
                REQUIRE(arrangement.cells.size() == 4);
                REQUIRE(count_num_halfedges(arrangement) == 15);
                REQUIRE(arrangement.vertices.size() == 9);
                REQUIRE(arrangement.unique_planes.size() == 6);
                validate_arrangement(arrangement);
            }

            SECTION("case 4")
            {
                planes.push_back({3, -1, -1});
                planes.push_back({-1, 3, -1});
                planes.push_back({-1, -1, 2});
                auto arrangement = compute_arrangement(planes);
                REQUIRE(arrangement.unique_plane_indices.size() == 6);
                REQUIRE(arrangement.cells.size() == 7);
                REQUIRE(count_num_halfedges(arrangement) == 27);
                REQUIRE(arrangement.vertices.size() == 12);
                REQUIRE(arrangement.unique_planes.size() == 6);
                validate_arrangement(arrangement);
            }

            SECTION("case 5")
            {
                planes.push_back({3, -1, -1});
                planes.push_back({-3, 1, 1});
                planes.push_back({-1, -1, 2});
                auto arrangement = compute_arrangement(planes);
                REQUIRE(arrangement.unique_plane_indices.size() == 6);
                REQUIRE(arrangement.cells.size() == 4);
                REQUIRE(count_num_halfedges(arrangement) == 15);
                REQUIRE(arrangement.vertices.size() == 8);
                REQUIRE(arrangement.unique_planes.size() == 5);
                validate_arrangement(arrangement);
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
        REQUIRE(arrangement.unique_planes.size() == 4);
        validate_arrangement(arrangement);
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
            REQUIRE(arrangement.unique_planes.size() == 5);
            validate_arrangement(arrangement);
        }
        SECTION("tri cross section")
        {
            planes.push_back({1, -1, -1, -1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 2);
            REQUIRE(count_num_half_faces(arrangement) == 9);
            REQUIRE(count_num_half_edges(arrangement) == 30);
            REQUIRE(arrangement.vertices.size() == 7);
            REQUIRE(arrangement.unique_planes.size() == 5);
            validate_arrangement(arrangement);
        }
        SECTION("cut through a vertex")
        {
            planes.push_back({0, -1, -1, 1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 2);
            REQUIRE(count_num_half_faces(arrangement) == 9);
            REQUIRE(count_num_half_edges(arrangement) == 28);
            REQUIRE(arrangement.vertices.size() == 6);
            REQUIRE(arrangement.unique_planes.size() == 5);
            validate_arrangement(arrangement);
        }
        SECTION("cut through an edge")
        {
            planes.push_back({0, 0, -1, 1});
            auto arrangement = compute_arrangement(planes);
            REQUIRE(count_num_cells(arrangement) == 2);
            REQUIRE(count_num_half_faces(arrangement) == 8);
            REQUIRE(count_num_half_edges(arrangement) == 24);
            REQUIRE(arrangement.vertices.size() == 5);
            REQUIRE(arrangement.unique_planes.size() == 5);
            validate_arrangement(arrangement);
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
            validate_arrangement(arrangement);
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
            REQUIRE(arrangement.unique_planes.size() == 6);
            validate_arrangement(arrangement);
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
            REQUIRE(arrangement.unique_planes.size() == 6);
            validate_arrangement(arrangement);
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
            REQUIRE(arrangement.unique_planes.size() == 6);
            validate_arrangement(arrangement);
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
            validate_arrangement(arrangement);
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
            REQUIRE(arrangement.unique_planes.size() == 6);
            validate_arrangement(arrangement);
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
            validate_arrangement(arrangement);
        }
    }
}

} // namespace

TEST_CASE("Arrangement 2D", "[arrangement][2D]")
{
    using namespace simplicial_arrangement;
    SECTION("Int") { test_2D<Int>(); }
    SECTION("double") { test_2D<double>(); }
}

TEST_CASE("Arrangement 3D", "[arrangement][3d]")
{
    using namespace simplicial_arrangement;
    SECTION("Int") { test_3D<Int>(); }
    SECTION("double") { test_3D<double>(); }
}
