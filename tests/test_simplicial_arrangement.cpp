#include <catch2/catch.hpp>

#include <simplicial_arrangement/SimplicialArrangement.h>

namespace {

template <typename Scalar>
void test_2D()
{
    using namespace simplicial_arrangement;

    std::function<size_t(const BSPNode<2>&)> count_num_cells;
    count_num_cells = [&](const BSPNode<2>& root) -> size_t {
        if (root.negative != nullptr && root.positive != nullptr) {
            return count_num_cells(*root.negative) + count_num_cells(*root.positive);
        } else {
            return 1;
        }
    };

    std::function<size_t(const BSPNode<2>&)> count_num_edges;
    count_num_edges = [&](const BSPNode<2>& root) -> size_t {
        if (root.negative != nullptr && root.positive != nullptr) {
            return count_num_edges(*root.negative) + count_num_edges(*root.positive);
        } else {
            return root.cell.edges.size();
        }
    };

    SimplicialArrangement<Scalar, 2> arrangement;
    std::vector<Plane<Scalar, 2>> planes;

    SECTION("0 implicit")
    {
        arrangement.initialize({});
        REQUIRE(arrangement.get_planes().size() == 3);
        REQUIRE(count_num_cells(arrangement.get_root()) == 1);
        REQUIRE(count_num_edges(arrangement.get_root()) == 3);
        REQUIRE(arrangement.get_vertex_count() == 3);
        REQUIRE(arrangement.get_num_unique_planes() == 3);
    }

    SECTION("1 implicit")
    {
        SECTION("normal case")
        {
            planes.push_back({1, -1, -1});
            arrangement.initialize(planes);
            REQUIRE(arrangement.get_planes().size() == 4);
            REQUIRE(count_num_cells(arrangement.get_root()) == 2);
            REQUIRE(count_num_edges(arrangement.get_root()) == 7);
            REQUIRE(arrangement.get_vertex_count() == 5);
            REQUIRE(arrangement.get_num_unique_planes() == 4);
        }
        SECTION("passing through one corner")
        {
            planes.push_back({1, -1, 0});
            arrangement.initialize(planes);
            REQUIRE(arrangement.get_planes().size() == 4);
            REQUIRE(count_num_cells(arrangement.get_root()) == 2);
            REQUIRE(count_num_edges(arrangement.get_root()) == 6);
            REQUIRE(arrangement.get_vertex_count() == 4);
            REQUIRE(arrangement.get_num_unique_planes() == 4);
        }
        SECTION("passing through two corners")
        {
            planes.push_back({1, 0, 0});
            arrangement.initialize(planes);
            REQUIRE(arrangement.get_planes().size() == 4);
            REQUIRE(count_num_cells(arrangement.get_root()) == 1);
            REQUIRE(count_num_edges(arrangement.get_root()) == 3);
            REQUIRE(arrangement.get_vertex_count() == 3);
            REQUIRE(arrangement.get_num_unique_planes() == 3);
        }
    }

    SECTION("2 implicits")
    {
        logger().set_level(spdlog::level::info);
        SECTION("not intersecting")
        {
            planes.push_back({1, -1, -1});
            planes.push_back({2, -1, -1});
            arrangement.initialize(planes);
            REQUIRE(count_num_cells(arrangement.get_root()) == 3);
            REQUIRE(count_num_edges(arrangement.get_root()) == 11);
            REQUIRE(arrangement.get_vertex_count() == 7);
            REQUIRE(arrangement.get_num_unique_planes() == 5);
        }

        SECTION("crossing")
        {
            planes.push_back({1, -1, -1});
            planes.push_back({-1, 2, -1});
            arrangement.initialize(planes);
            REQUIRE(count_num_cells(arrangement.get_root()) == 4);
            REQUIRE(count_num_edges(arrangement.get_root()) == 15);
            REQUIRE(arrangement.get_vertex_count() == 8);
            REQUIRE(arrangement.get_num_unique_planes() == 5);
        }

        SECTION("touching")
        {
            planes.push_back({1, -1, -1});
            planes.push_back({-1, 1, -1});
            arrangement.initialize(planes);
            REQUIRE(count_num_cells(arrangement.get_root()) == 3);
            REQUIRE(count_num_edges(arrangement.get_root()) == 10);
            REQUIRE(arrangement.get_vertex_count() == 6);
            REQUIRE(arrangement.get_num_unique_planes() == 5);
        }

        SECTION("overlapping")
        {
            planes.push_back({1, -1, -1});
            planes.push_back({1, -1, -1});
            arrangement.initialize(planes);
            REQUIRE(count_num_cells(arrangement.get_root()) == 2);
            REQUIRE(count_num_edges(arrangement.get_root()) == 7);
            REQUIRE(arrangement.get_vertex_count() == 5);
            REQUIRE(arrangement.get_num_unique_planes() == 4);
        }
    }
    SECTION("3 implicits")
    {
        SECTION("non-intersecting")
        {
            planes.push_back({1, -1, -1});
            planes.push_back({2, -1, -1});
            planes.push_back({3, -1, -1});
            arrangement.initialize(planes);
            REQUIRE(count_num_cells(arrangement.get_root()) == 4);
            REQUIRE(count_num_edges(arrangement.get_root()) == 15);
            REQUIRE(arrangement.get_vertex_count() == 9);
            REQUIRE(arrangement.get_num_unique_planes() == 6);
        }

        SECTION("intersecting")
        {
            SECTION("case 0")
            {
                planes.push_back({3, -1, -1});
                planes.push_back({-1, 3, -1});
                planes.push_back({-1, -1, 3});
                arrangement.initialize(planes);
                REQUIRE(count_num_cells(arrangement.get_root()) == 7);
                REQUIRE(count_num_edges(arrangement.get_root()) == 27);
                REQUIRE(arrangement.get_vertex_count() == 12);
                REQUIRE(arrangement.get_num_unique_planes() == 6);
            }

            SECTION("case 1")
            {
                planes.push_back({2, -1, -1});
                planes.push_back({-1, 2, -1});
                planes.push_back({-1, -1, 2});
                arrangement.initialize(planes);
                REQUIRE(count_num_cells(arrangement.get_root()) == 6);
                REQUIRE(count_num_edges(arrangement.get_root()) == 21);
                REQUIRE(arrangement.get_vertex_count() == 10);
                REQUIRE(arrangement.get_num_unique_planes() == 6);
            }

            SECTION("case 2")
            {
                planes.push_back({1, -1, -1});
                planes.push_back({-1, 1, -1});
                planes.push_back({-1, -1, 1});
                arrangement.initialize(planes);
                REQUIRE(count_num_cells(arrangement.get_root()) == 4);
                REQUIRE(count_num_edges(arrangement.get_root()) == 12);
                REQUIRE(arrangement.get_vertex_count() == 6);
                REQUIRE(arrangement.get_num_unique_planes() == 6);
            }

            SECTION("case 3")
            {
                planes.push_back({1, -2, -2});
                planes.push_back({-2, 1, -2});
                planes.push_back({-2, -2, 1});
                arrangement.initialize(planes);
                REQUIRE(count_num_cells(arrangement.get_root()) == 4);
                REQUIRE(count_num_edges(arrangement.get_root()) == 15);
                REQUIRE(arrangement.get_vertex_count() == 9);
                REQUIRE(arrangement.get_num_unique_planes() == 6);
            }

            SECTION("case 4")
            {
                planes.push_back({3, -1, -1});
                planes.push_back({-1, 3, -1});
                planes.push_back({-1, -1, 2});
                arrangement.initialize(planes);
                REQUIRE(count_num_cells(arrangement.get_root()) == 7);
                REQUIRE(count_num_edges(arrangement.get_root()) == 27);
                REQUIRE(arrangement.get_vertex_count() == 12);
                REQUIRE(arrangement.get_num_unique_planes() == 6);
            }

            SECTION("case 5")
            {
                planes.push_back({3, -1, -1});
                planes.push_back({-3, 1, 1});
                planes.push_back({-1, -1, 2});
                arrangement.initialize(planes);
                REQUIRE(count_num_cells(arrangement.get_root()) == 4);
                REQUIRE(count_num_edges(arrangement.get_root()) == 15);
                REQUIRE(arrangement.get_vertex_count() == 8);
                REQUIRE(arrangement.get_num_unique_planes() == 5);
            }
        }
    }
}

template <typename Scalar>
void test_3D()
{
    using namespace simplicial_arrangement;

    std::function<size_t(const BSPNode<3>&)> count_num_cells;
    count_num_cells = [&](const BSPNode<3>& root) -> size_t {
        if (root.negative != nullptr && root.positive != nullptr) {
            return count_num_cells(*root.negative) + count_num_cells(*root.positive);
        } else {
            return 1;
        }
    };

    std::function<size_t(const BSPNode<3>&)> count_num_faces;
    count_num_faces = [&](const BSPNode<3>& root) -> size_t {
        if (root.negative != nullptr && root.positive != nullptr) {
            return count_num_faces(*root.negative) + count_num_faces(*root.positive);
        } else {
            return root.cell.faces.size();
        }
    };

    std::function<size_t(const BSPNode<3>&)> count_num_half_edges;
    count_num_half_edges = [&](const BSPNode<3>& root) -> size_t {
        if (root.negative != nullptr && root.positive != nullptr) {
            return count_num_half_edges(*root.negative) + count_num_half_edges(*root.positive);
        } else {
            size_t num_half_edges = 0;
            std::for_each(root.cell.faces.begin(), root.cell.faces.end(), [&](const auto& face) {
                num_half_edges += face.edge_planes.size();
            });
            return num_half_edges;
        }
    };

    SimplicialArrangement<Scalar, 3> arrangement;
    std::vector<Plane<Scalar, 3>> planes;

    SECTION("0 implicit")
    {
        arrangement.initialize(planes);
        REQUIRE(count_num_cells(arrangement.get_root()) == 1);
        REQUIRE(count_num_faces(arrangement.get_root()) == 4);
        REQUIRE(count_num_half_edges(arrangement.get_root()) == 12);
        REQUIRE(arrangement.get_vertex_count() == 4);
        REQUIRE(arrangement.get_num_unique_planes() == 4);
    }

    SECTION("1 implicit")
    {
        SECTION("quad cross section")
        {
            planes.push_back({1, 1, -1, -1});
            arrangement.initialize(planes);
            REQUIRE(count_num_cells(arrangement.get_root()) == 2);
            REQUIRE(count_num_faces(arrangement.get_root()) == 10);
            REQUIRE(count_num_half_edges(arrangement.get_root()) == 36);
            REQUIRE(arrangement.get_vertex_count() == 8);
            REQUIRE(arrangement.get_num_unique_planes() == 5);
        }
        SECTION("tri cross section")
        {
            planes.push_back({1, -1, -1, -1});
            arrangement.initialize(planes);
            REQUIRE(count_num_cells(arrangement.get_root()) == 2);
            REQUIRE(count_num_faces(arrangement.get_root()) == 9);
            REQUIRE(count_num_half_edges(arrangement.get_root()) == 30);
            REQUIRE(arrangement.get_vertex_count() == 7);
            REQUIRE(arrangement.get_num_unique_planes() == 5);
        }
        SECTION("cut through a vertex")
        {
            planes.push_back({0, -1, -1, 1});
            arrangement.initialize(planes);
            REQUIRE(count_num_cells(arrangement.get_root()) == 2);
            REQUIRE(count_num_faces(arrangement.get_root()) == 9);
            REQUIRE(count_num_half_edges(arrangement.get_root()) == 28);
            REQUIRE(arrangement.get_vertex_count() == 6);
            REQUIRE(arrangement.get_num_unique_planes() == 5);
        }
        SECTION("cut through an edge")
        {
            planes.push_back({0, 0, -1, 1});
            arrangement.initialize(planes);
            REQUIRE(count_num_cells(arrangement.get_root()) == 2);
            REQUIRE(count_num_faces(arrangement.get_root()) == 8);
            REQUIRE(count_num_half_edges(arrangement.get_root()) == 24);
            REQUIRE(arrangement.get_vertex_count() == 5);
            REQUIRE(arrangement.get_num_unique_planes() == 5);
        }
        SECTION("cut through a face")
        {
            planes.push_back({0, 0, 0, 1});
            arrangement.initialize(planes);
            REQUIRE(count_num_cells(arrangement.get_root()) == 1);
            REQUIRE(count_num_faces(arrangement.get_root()) == 4);
            REQUIRE(count_num_half_edges(arrangement.get_root()) == 12);
            REQUIRE(arrangement.get_vertex_count() == 4);
            REQUIRE(arrangement.get_num_unique_planes() == 4);
        }
    }

    SECTION("2 implicits")
    {
        SECTION("not intersecting cuts")
        {
            planes.push_back({1, -2, -2, -2});
            planes.push_back({-2, 1, -2, -2});
            arrangement.initialize(planes);
            REQUIRE(count_num_cells(arrangement.get_root()) == 3);
            REQUIRE(count_num_faces(arrangement.get_root()) == 14);
            REQUIRE(count_num_half_edges(arrangement.get_root()) == 48);
            REQUIRE(arrangement.get_vertex_count() == 10);
            REQUIRE(arrangement.get_num_unique_planes() == 6);
        }
        SECTION("crossing")
        {
            planes.push_back({2, -1, -1, -1});
            planes.push_back({-1, 2, -1, -1});
            arrangement.initialize(planes);
            REQUIRE(count_num_cells(arrangement.get_root()) == 4);
            REQUIRE(count_num_faces(arrangement.get_root()) == 20);
            REQUIRE(count_num_half_edges(arrangement.get_root()) == 72);
            REQUIRE(arrangement.get_vertex_count() == 12);
            REQUIRE(arrangement.get_num_unique_planes() == 6);
        }
        SECTION("touching")
        {
            planes.push_back({1, -1, -1, -1});
            planes.push_back({-1, 1, -1, -1});
            arrangement.initialize(planes);
            REQUIRE(count_num_cells(arrangement.get_root()) == 3);
            REQUIRE(count_num_faces(arrangement.get_root()) == 14);
            REQUIRE(count_num_half_edges(arrangement.get_root()) == 46);
            REQUIRE(arrangement.get_vertex_count() == 9);
            REQUIRE(arrangement.get_num_unique_planes() == 6);
        }
        SECTION("duplicated cuts")
        {
            planes.push_back({1, 1, -1, -1});
            planes.push_back({2, 2, -2, -2});
            arrangement.initialize(planes);
            REQUIRE(count_num_cells(arrangement.get_root()) == 2);
            REQUIRE(count_num_faces(arrangement.get_root()) == 10);
            REQUIRE(count_num_half_edges(arrangement.get_root()) == 36);
            REQUIRE(arrangement.get_vertex_count() == 8);
            REQUIRE(arrangement.get_num_unique_planes() == 5);
        }
        SECTION("2 cuts through the same vertex")
        {
            planes.push_back({0, 0, 1, -1});
            planes.push_back({0, 2, -2, 0});
            arrangement.initialize(planes);
            REQUIRE(count_num_cells(arrangement.get_root()) == 4);
            REQUIRE(count_num_faces(arrangement.get_root()) == 17);
            REQUIRE(count_num_half_edges(arrangement.get_root()) == 52);
            REQUIRE(arrangement.get_vertex_count() == 7);
            REQUIRE(arrangement.get_num_unique_planes() == 6);
        }
        SECTION("2 cuts through faces")
        {
            planes.push_back({0, 0, 0, 1});
            planes.push_back({0, 0, 0, 2});
            arrangement.initialize(planes);
            REQUIRE(count_num_cells(arrangement.get_root()) == 1);
            REQUIRE(count_num_faces(arrangement.get_root()) == 4);
            REQUIRE(count_num_half_edges(arrangement.get_root()) == 12);
            REQUIRE(arrangement.get_vertex_count() == 4);
            REQUIRE(arrangement.get_num_unique_planes() == 4);
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
