#include <simplicial_arrangement/material_interface.h>
#include <catch2/catch.hpp>

#include <implicit_predicates/implicit_predicates.h>

namespace {
template <typename Scalar>
void test_2D()
{
    using namespace simplicial_arrangement;

    std::vector<Material<Scalar, 2>> materials;

    SECTION("1 implicit")
    {
        materials.push_back({1, 1, 1});
        auto mi = compute_material_interface(materials);
        REQUIRE(mi.cells.size() == 1);
        REQUIRE(mi.faces.size() == 3);
        REQUIRE(mi.vertices.size() == 3);
    }

    SECTION("2 implicits")
    {
        materials.push_back({1, 0, 0});
        materials.push_back({0, 1, 0});
        auto mi = compute_material_interface(materials);
        REQUIRE(mi.cells.size() == 2);
        REQUIRE(mi.faces.size() == 5);
        REQUIRE(mi.vertices.size() == 4);
    }
}
} // namespace

TEST_CASE("Material interface 2D", "[material_interface][2D]")
{
    using namespace simplicial_arrangement;
    SECTION("Int") { test_2D<Int>(); }
    SECTION("double") { test_2D<double>(); }
}
