#include <simplicial_arrangement/lookup_table.h>
#include <simplicial_arrangement/simplicial_arrangement.h>

#include <lookup_table_forward_declarations.h>

#include <catch2/catch.hpp>
#include <spdlog/spdlog.h>

TEST_CASE("Lookup table", "[lookup]")
{
    using namespace simplicial_arrangement;

    bool loaded = load_lookup_table();
    REQUIRE(loaded);

    REQUIRE(one_func_lookup_table->size() == 16);
    REQUIRE(two_func_lookup_table->size() == 256);
}

TEST_CASE("Material Interface lookup", "[mi][lookup]")
{
    using namespace simplicial_arrangement;
    std::vector<Material<Int, 3>> materials;
    REQUIRE(load_lookup_table());
    spdlog::set_level(spdlog::level::info);

    auto assert_same = [](const MaterialInterface<3>& mi0, const MaterialInterface<3>& mi1) {
        REQUIRE(mi0.vertices.size() == mi1.vertices.size());
        REQUIRE(mi0.faces.size() == mi1.faces.size());
        REQUIRE(mi0.cells.size() == mi1.cells.size());
    };

    SECTION("Case 1")
    {
        materials.push_back({1, 1, 1, 1});
        materials.push_back({2, 2, 2, 2});
        materials.push_back({3, 3, 3, 3});

        enable_lookup_table();
        const auto mi0 = compute_material_interface(materials);
        disable_lookup_table();
        const auto mi1 = compute_material_interface(materials);

        assert_same(mi0, mi1);
    }
    SECTION("Case 2")
    {
        materials.push_back({2, 1, 1, 1});
        materials.push_back({1, 2, 2, 2});
        materials.push_back({3, 3, 3, 3});

        enable_lookup_table();
        const auto mi0 = compute_material_interface(materials);
        disable_lookup_table();
        const auto mi1 = compute_material_interface(materials);

        assert_same(mi0, mi1);
    }
    SECTION("Case 3")
    {
        materials.push_back({1, 1, 1, 1});
        materials.push_back({2, 2, 2, 2});
        materials.push_back({0, 3, 3, 3});

        enable_lookup_table();
        const auto mi0 = compute_material_interface(materials);
        disable_lookup_table();
        const auto mi1 = compute_material_interface(materials);

        assert_same(mi0, mi1);
    }
    SECTION("Case 4")
    {
        materials.push_back({3, 3, 1, 1});
        materials.push_back({1, 1, 2, 2});
        materials.push_back({0, 0, 0, 0});

        enable_lookup_table();
        const auto mi0 = compute_material_interface(materials);
        disable_lookup_table();
        const auto mi1 = compute_material_interface(materials);

        assert_same(mi0, mi1);
    }
    SECTION("Case 5")
    {
        materials.push_back({1, 1, 1, 2});
        materials.push_back({2, 2, 2, 1});
        materials.push_back({0, 0, 0, 0});

        enable_lookup_table();
        const auto mi0 = compute_material_interface(materials);
        disable_lookup_table();
        const auto mi1 = compute_material_interface(materials);

        assert_same(mi0, mi1);
    }
    SECTION("Case 6")
    {
        materials.push_back({1, 4, 1, 2});
        materials.push_back({2, 2, 2, 1});
        materials.push_back({4, 1, 0, 0});

        enable_lookup_table();
        const auto mi0 = compute_material_interface(materials);
        disable_lookup_table();
        const auto mi1 = compute_material_interface(materials);

        assert_same(mi0, mi1);
    }
    SECTION("Case 7")
    {
        materials.push_back({1, 4, 3, 1});
        materials.push_back({2, 2, 2, 2});
        materials.push_back({4, 1, 0, 4});

        enable_lookup_table();
        const auto mi0 = compute_material_interface(materials);
        disable_lookup_table();
        const auto mi1 = compute_material_interface(materials);

        assert_same(mi0, mi1);
    }
}
