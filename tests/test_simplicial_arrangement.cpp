#include <catch2/catch.hpp>

#include <simplicial_arrangement/SimplicialArrangement.h>

TEST_CASE("Arrangement 2D", "[arrangement][2D]")
{
    using namespace simplicial_arrangement;
    using Scalar = Int;

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
    }

    SECTION("1 implicit")
    {
        SECTION("normal case") {
            planes.push_back({1, -1, -1});
            arrangement.initialize(planes);
            REQUIRE(arrangement.get_planes().size() == 4);
            REQUIRE(count_num_cells(arrangement.get_root()) == 2);
            REQUIRE(count_num_edges(arrangement.get_root()) == 7);
        }
        SECTION("passing through one corner") {
            planes.push_back({1, -1, 0});
            arrangement.initialize(planes);
            REQUIRE(arrangement.get_planes().size() == 4);
            REQUIRE(count_num_cells(arrangement.get_root()) == 2);
            REQUIRE(count_num_edges(arrangement.get_root()) == 6);
        }
        SECTION("passing through two corners") {
            planes.push_back({1, 0, 0});
            arrangement.initialize(planes);
            REQUIRE(arrangement.get_planes().size() == 4);
            REQUIRE(count_num_cells(arrangement.get_root()) == 1);
            REQUIRE(count_num_edges(arrangement.get_root()) == 3);
        }
    }

    SECTION("2 implicits")
    {
        SECTION("not intersecting")
        {
            planes.push_back({1, -1, -1});
            planes.push_back({2, -1, -1});
            arrangement.initialize(planes);
            REQUIRE(count_num_cells(arrangement.get_root()) == 3);
            REQUIRE(count_num_edges(arrangement.get_root()) == 11);
        }

        SECTION("crossing")
        {
            planes.push_back({1, -1, -1});
            planes.push_back({-1, 2, -1});
            arrangement.initialize(planes);
            REQUIRE(count_num_cells(arrangement.get_root()) == 4);
            REQUIRE(count_num_edges(arrangement.get_root()) == 15);
        }

        SECTION("touching")
        {
            planes.push_back({1, -1, -1});
            planes.push_back({-1, 1, -1});
            arrangement.initialize(planes);
            REQUIRE(count_num_cells(arrangement.get_root()) == 3);
            REQUIRE(count_num_edges(arrangement.get_root()) == 10);
        }

        SECTION("overlapping")
        {
            planes.push_back({1, -1, -1});
            planes.push_back({1, -1, -1});
            arrangement.initialize(planes);
            REQUIRE(count_num_cells(arrangement.get_root()) == 2);
            REQUIRE(count_num_edges(arrangement.get_root()) == 7);
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
            }

            SECTION("case 1")
            {
                planes.push_back({2, -1, -1});
                planes.push_back({-1, 2, -1});
                planes.push_back({-1, -1, 2});
                arrangement.initialize(planes);
                REQUIRE(count_num_cells(arrangement.get_root()) == 6);
                REQUIRE(count_num_edges(arrangement.get_root()) == 21);
            }

            SECTION("case 2")
            {
                planes.push_back({1, -1, -1});
                planes.push_back({-1, 1, -1});
                planes.push_back({-1, -1, 1});
                arrangement.initialize(planes);
                REQUIRE(count_num_cells(arrangement.get_root()) == 4);
                REQUIRE(count_num_edges(arrangement.get_root()) == 12);
            }

            SECTION("case 3")
            {
                planes.push_back({1, -2, -2});
                planes.push_back({-2, 1, -2});
                planes.push_back({-2, -2, 1});
                arrangement.initialize(planes);
                REQUIRE(count_num_cells(arrangement.get_root()) == 4);
                REQUIRE(count_num_edges(arrangement.get_root()) == 15);
            }

            SECTION("case 4")
            {
                planes.push_back({3, -1, -1});
                planes.push_back({-1, 3, -1});
                planes.push_back({-1, -1, 2});
                arrangement.initialize(planes);
                REQUIRE(count_num_cells(arrangement.get_root()) == 7);
                REQUIRE(count_num_edges(arrangement.get_root()) == 27);
            }
        }
    }
}
