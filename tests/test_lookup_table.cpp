#include <simplicial_arrangement/lookup_table.h>
#include <simplicial_arrangement/simplicial_arrangement.h>

#include <lookup_table_forward_declarations.h>

#include <absl/container/flat_hash_map.h>

#include <spdlog/spdlog.h>
#include <catch2/catch.hpp>

TEST_CASE("Simplicial arrangement lookup", "[ar][lookup]")
{
    using namespace simplicial_arrangement;
    std::vector<Plane<double, 3>> planes;
    REQUIRE(load_lookup_table(ARRANGEMENT));
    spdlog::set_level(spdlog::level::info);

    auto assert_same = [](const Arrangement<3>& ar0, const Arrangement<3>& ar1) {
        REQUIRE(ar0.vertices.size() == ar1.vertices.size());
        REQUIRE(ar0.faces.size() == ar1.faces.size());
        REQUIRE(ar0.cells.size() == ar1.cells.size());

        auto count_vertices_on_bd_face = [](const auto& ar) {
            std::array<size_t, 4> count{0, 0, 0, 0};
            for (const auto& v : ar.vertices) {
                for (size_t i = 0; i < 4; i++) {
                    if (v[0] == i || v[1] == i || v[2] == i) {
                        count[i]++;
                    }
                }
            }
            return count;
        };

        auto count0 = count_vertices_on_bd_face(ar0);
        auto count1 = count_vertices_on_bd_face(ar1);
        REQUIRE(count0 == count1);

        absl::flat_hash_map<std::pair<size_t, size_t>, size_t> cell_to_plane;

        for (const auto& face : ar0.faces) {
            if (face.supporting_plane < 4) continue; // Boundary face.
            std::pair<size_t, size_t> key{face.positive_cell, face.negative_cell};
            cell_to_plane.emplace(key, face.supporting_plane);
        }

        for (const auto& face : ar1.faces) {
            if (face.supporting_plane < 4) continue; // Boundary face.
            std::pair<size_t, size_t> key{face.positive_cell, face.negative_cell};
            const auto itr = cell_to_plane.find(key);
            REQUIRE(itr != cell_to_plane.end());
            REQUIRE(itr->second == face.supporting_plane);
        }
    };

    SECTION("Case 1")
    {
        planes.push_back({1, 1, 1, 1});
        planes.push_back({0, 0, 0, 0});
    }
    SECTION("Case 2")
    {
        planes.push_back({1, -1, -1, -1});
        planes.push_back({-1, 1, -1, -1});
    }
    SECTION("Case 3")
    {
        planes.push_back({1, -1, -1, -1});
        planes.push_back({2, -1, -1, -1});
    }
    SECTION("Case 4")
    {
        planes.push_back({1, -1, -1, -1});
        planes.push_back({-1, 2, -1, -1});
    }
    SECTION("Case 5")
    {
        planes.push_back({1, 1, -1, -1});
        planes.push_back({-1, 2, -1, -1});
    }
    SECTION("Case 6")
    {
        planes.push_back({1, -1, -1, -1});
        planes.push_back({2, -2, -2, -2});
    }
    SECTION("Case 7")
    {
        planes.push_back({2, -1, -1, -1});
        planes.push_back({-2, 4, -2, -2});
    }
    SECTION("Case 8")
    {
        planes.push_back({2, 2, -1, -1});
        planes.push_back({-2, 5, 6, -2});
    }
    SECTION("Case 9") { planes.push_back({1, 1, 1, -1}); }
    SECTION("Case 10") { planes.push_back({1, 1, -2, -1}); }
    SECTION("Case 11")
    {
        planes.push_back({-0.011580921511964659,
            -0.093672652278739665,
            0.025375829081217716,
            0.0066544614307747496});
        planes.push_back({-0.021899454112778494,
            -0.017737799597463133,
            0.014110607076315729,
            0.10106216401814538});
    }

    enable_lookup_table();
    const auto ar0 = compute_arrangement(planes);
    disable_lookup_table();
    const auto ar1 = compute_arrangement(planes);

    assert_same(ar0, ar1);
}

TEST_CASE("Material Interface lookup", "[mi][lookup]")
{
    using namespace simplicial_arrangement;
    std::vector<Material<double, 3>> materials;
    REQUIRE(load_lookup_table(MATERIAL_INTERFACE));
    spdlog::set_level(spdlog::level::info);

    auto assert_same = [](const MaterialInterface<3>& mi0, const MaterialInterface<3>& mi1) {
        REQUIRE(mi0.vertices.size() == mi1.vertices.size());
        REQUIRE(mi0.faces.size() == mi1.faces.size());
        REQUIRE(mi0.cells.size() == mi1.cells.size());

        auto count_vertices_on_bd_face = [](const auto& mi) {
            std::array<size_t, 4> count{0, 0, 0, 0};
            for (const auto& v : mi.vertices) {
                for (size_t i = 0; i < 4; i++) {
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
        for (size_t i = 0; i < num_cells; i++) {
            const auto& c = mi0.cells[i];
            material_to_cell_map.insert({c.material_label, i});
        }
        for (size_t i = 0; i < num_cells; i++) {
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
    }
    SECTION("Case 2")
    {
        materials.push_back({2, 1, 1, 1});
        materials.push_back({1, 2, 2, 2});
        materials.push_back({3, 3, 3, 3});
    }
    SECTION("Case 3")
    {
        materials.push_back({1, 1, 1, 1});
        materials.push_back({2, 2, 2, 2});
        materials.push_back({0, 3, 3, 3});
    }
    SECTION("Case 4")
    {
        materials.push_back({3, 3, 1, 1});
        materials.push_back({1, 1, 2, 2});
        materials.push_back({0, 0, 0, 0});
    }
    SECTION("Case 5")
    {
        materials.push_back({1, 1, 1, 2});
        materials.push_back({2, 2, 2, 1});
        materials.push_back({0, 0, 0, 0});
    }
    SECTION("Case 6")
    {
        materials.push_back({1, 4, 1, 2});
        materials.push_back({2, 2, 2, 1});
        materials.push_back({4, 1, 0, 0});
    }
    SECTION("Case 7")
    {
        materials.push_back({1, 4, 3, 1});
        materials.push_back({2, 2, 2, 2});
        materials.push_back({4, 1, 0, 4});
    }

    enable_lookup_table();
    const auto mi0 = compute_material_interface(materials);
    disable_lookup_table();
    const auto mi1 = compute_material_interface(materials);

    assert_same(mi0, mi1);
}
