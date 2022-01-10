#include <simplicial_arrangement/lookup_table.h>
#include <simplicial_arrangement/simplicial_arrangement.h>

#include <lookup_table_forward_declarations.h>

#include <absl/container/flat_hash_map.h>

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

        auto count_vertices_on_bd_face = [](const auto& mi) {
            std::array<size_t, 4> count{0, 0, 0, 0};
            for (const auto& v : mi.vertices) {
                for (size_t i=0; i<4; i++) {
                    if (v[0] == i || v[1] == i || v[2] == i || v[3] == i) {
                        count[i]++;
                    }
                }
            }
            return count;
        };

        std::array<size_t, 4> count0 = count_vertices_on_bd_face(mi0);
        std::array<size_t, 4> count1 = count_vertices_on_bd_face(mi1);
        REQUIRE(count0 == count1);

        absl::flat_hash_map<size_t, size_t> material_to_cell_map;

        const size_t num_cells = mi0.cells.size();
        for (size_t i=0; i<num_cells; i++) {
            const auto& c = mi0.cells[i];
            material_to_cell_map.insert({c.material_label, i});
        }
        for (size_t i=0; i<num_cells; i++) {
            const auto& c1 = mi1.cells[i];
            const size_t j = material_to_cell_map[c1.material_label];
            const auto& c0 = mi0.cells[j];

            REQUIRE(c0.faces.size() == c1.faces.size());
        }
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
