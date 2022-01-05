#include <simplicial_arrangement/material_interface.h>
#include <spdlog/spdlog.h>
#include <catch2/catch.hpp>

#include <implicit_predicates/implicit_predicates.h>

namespace {

template<int DIM>
void validate(const simplicial_arrangement::MaterialInterface<DIM>& mi)
{
    const size_t num_vertices = mi.vertices.size();
    const size_t num_faces = mi.faces.size();
    const size_t num_cells = mi.cells.size();
    REQUIRE(num_cells > 0);

    std::vector<size_t> vertex_valence(num_vertices, 0);
    std::vector<size_t> face_valence(num_faces, 0);

    for (const auto& f : mi.faces) {
        for (const auto vi : f.vertices) {
            REQUIRE(vi < num_vertices);
            vertex_valence[vi]++;
        }

        REQUIRE(f.positive_material_label != f.negative_material_label);
    }
    for (const auto& c : mi.cells) {
        for (const auto fi : c.faces) {
            REQUIRE(fi < num_faces);
            face_valence[fi]++;

            const auto& f = mi.faces[fi];
            REQUIRE(((c.material_label == f.positive_material_label) ||
                     (c.material_label == f.negative_material_label)));
        }
    }
    REQUIRE(
        std::all_of(vertex_valence.begin(), vertex_valence.end(), [](size_t v) { return v > 0; }));
    REQUIRE(std::all_of(face_valence.begin(), face_valence.end(), [](size_t v) { return v > 0; }));
}

template <typename Scalar>
void test_2D()
{
    using namespace simplicial_arrangement;

    std::vector<Material<Scalar, 2>> materials;
    spdlog::set_level(spdlog::level::info);

    SECTION("1 implicit")
    {
        materials.push_back({1, 1, 1});
        auto mi = compute_material_interface(materials);
        REQUIRE(mi.cells.size() == 1);
        REQUIRE(mi.faces.size() == 3);
        REQUIRE(mi.vertices.size() == 3);
        validate(mi);
    }

    SECTION("2 implicits")
    {
        SECTION("Case 1")
        {
            materials.push_back({1, 0, 0});
            materials.push_back({0, 1, 0});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 2);
            REQUIRE(mi.faces.size() == 5);
            REQUIRE(mi.vertices.size() == 4);
            validate(mi);
        }
        SECTION("Case 2")
        {
            materials.push_back({1, 0, 0});
            materials.push_back({0, 2, 1});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 2);
            REQUIRE(mi.faces.size() == 6);
            REQUIRE(mi.vertices.size() == 5);
            validate(mi);
        }
    }

    SECTION("3 implicits")
    {
        SECTION("Case 1")
        {
            materials.push_back({1, 0, 0});
            materials.push_back({0, 1, 0});
            materials.push_back({0, 0, 1});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 3);
            REQUIRE(mi.faces.size() == 9);
            REQUIRE(mi.vertices.size() == 7);
            validate(mi);
        }
    }

    SECTION("4 implicits")
    {
        SECTION("Case 1")
        {
            materials.push_back({6, 0, 0});
            materials.push_back({0, 6, 0});
            materials.push_back({0, 0, 6});
            materials.push_back({3, 3, 3});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 4);
            REQUIRE(mi.faces.size() == 9);
            REQUIRE(mi.vertices.size() == 6);
            validate(mi);
        }
        SECTION("Case 2")
        {
            materials.push_back({12, 0, 0});
            materials.push_back({0, 12, 0});
            materials.push_back({0, 0, 12});
            materials.push_back({5, 5, 5});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 4);
            REQUIRE(mi.faces.size() == 12);
            REQUIRE(mi.vertices.size() == 9);
            validate(mi);
        }
        SECTION("Case 3: just toucing")
        {
            materials.push_back({12, 0, 0});
            materials.push_back({0, 12, 0});
            materials.push_back({0, 0, 12});
            materials.push_back({4, 4, 4});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 3);
            REQUIRE(mi.faces.size() == 9);
            REQUIRE(mi.vertices.size() == 7);
            validate(mi);
        }
        SECTION("Case 4")
        {
            materials.push_back({12, 0, 0});
            materials.push_back({0, 12, 0});
            materials.push_back({0, 0, 12});
            materials.push_back({10, 10, 10});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 4);
            REQUIRE(mi.faces.size() == 12);
            REQUIRE(mi.vertices.size() == 9);
            validate(mi);
        }
    }
}

template <typename Scalar>
void test_3D()
{
    using namespace simplicial_arrangement;

    std::vector<Material<Scalar, 3>> materials;
    spdlog::set_level(spdlog::level::info);

    SECTION("1 implicit")
    {
        materials.push_back({1, 1, 1, 1});
        auto mi = compute_material_interface(materials);
        REQUIRE(mi.cells.size() == 1);
        REQUIRE(mi.faces.size() == 4);
        REQUIRE(mi.vertices.size() == 4);
        validate(mi);
    }

    SECTION("2 implicits")
    {
        SECTION("Case 1")
        {
            materials.push_back({1, 0, 0, 0});
            materials.push_back({0, 1, 0, 0});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 2);
            REQUIRE(mi.faces.size() == 7);
            REQUIRE(mi.vertices.size() == 5);
            validate(mi);
        }
        SECTION("Case 2")
        {
            materials.push_back({1, 0, 0, 0});
            materials.push_back({-1, 1, -1, -1});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 2);
            REQUIRE(mi.faces.size() == 8);
            REQUIRE(mi.vertices.size() == 7);
            validate(mi);
        }
        SECTION("Case 3")
        {
            materials.push_back({1, 1, 0, 0});
            materials.push_back({0, 0, 1, 1});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 2);
            REQUIRE(mi.faces.size() == 9);
            REQUIRE(mi.vertices.size() == 8);
            validate(mi);
        }
        SECTION("Case 4")
        {
            materials.push_back({1, 0, 0, 0});
            materials.push_back({0, 0, 0, 0});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 1);
            REQUIRE(mi.faces.size() == 4);
            REQUIRE(mi.vertices.size() == 4);
            validate(mi);
        }
        SECTION("Case 5")
        {
            materials.push_back({1, 0, 0, 0});
            materials.push_back({1, 1, 1, -1});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 2);
            REQUIRE(mi.faces.size() == 8);
            REQUIRE(mi.vertices.size() == 6);
            validate(mi);
        }
        SECTION("Case 6")
        {
            materials.push_back({1, 1, 0, 0});
            materials.push_back({0, 0, 1, 1});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 2);
            REQUIRE(mi.faces.size() == 9);
            REQUIRE(mi.vertices.size() == 8);
            validate(mi);
        }
    }
    SECTION("3 implicits") {
        SECTION("Case 1") {
            materials.push_back({1, 0, 0, 0});
            materials.push_back({0, 1, 0, 0});
            materials.push_back({0, 0, 1, 0});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 3);
            REQUIRE(mi.faces.size() == 12);
            REQUIRE(mi.vertices.size() == 8);
            validate(mi);
        }
        SECTION("Case 2") {
            materials.push_back({1, 0, 0, 0});
            materials.push_back({0, 1, 0, 0});
            materials.push_back({2, 0, 0, 0});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 2);
            REQUIRE(mi.faces.size() == 7);
            REQUIRE(mi.vertices.size() == 5);
            validate(mi);
        }
        SECTION("Case 3") {
            materials.push_back({1, 0, 0, 0});
            materials.push_back({0, 1, 0, 0});
            materials.push_back({0, 0, 1, 1});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 3);
            REQUIRE(mi.faces.size() == 13);
            REQUIRE(mi.vertices.size() == 11);
            validate(mi);
        }
    }
    SECTION("4 implicits") {
        SECTION("Case 1") {
            materials.push_back({1, 0, 0, 0});
            materials.push_back({0, 1, 0, 0});
            materials.push_back({0, 0, 1, 0});
            materials.push_back({0, 0, 0, 1});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 4);
            REQUIRE(mi.faces.size() == 18);
            REQUIRE(mi.vertices.size() == 15);
            validate(mi);
        }
        SECTION("Case 2") {
            materials.push_back({1, 0, 0, 0});
            materials.push_back({0, 1, 0, 0});
            materials.push_back({0, 0, 1, 0});
            materials.push_back({5, 0, 0, 0});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 3);
            REQUIRE(mi.faces.size() == 12);
            REQUIRE(mi.vertices.size() == 8);
            validate(mi);
        }
        SECTION("Case 3") {
            materials.push_back({5, 0, 0, 0});
            materials.push_back({0, 1, 0, 0});
            materials.push_back({0, 0, 1, 0});
            materials.push_back({-5, 0, 0, 0});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 3);
            REQUIRE(mi.faces.size() == 12);
            REQUIRE(mi.vertices.size() == 8);
            validate(mi);
        }
        SECTION("Case 4") {
            materials.push_back({1, 1, 0, 0});
            materials.push_back({0, 1, 1, 0});
            materials.push_back({0, 0, 1, 1});
            materials.push_back({1, 0, 0, 1});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 4);
            REQUIRE(mi.faces.size() == 12);
            REQUIRE(mi.vertices.size() == 6);
            validate(mi);
        }
        SECTION("Case 5") {
            // One material above all others.
            materials.push_back({1, 1, 0, 0});
            materials.push_back({0, 1, 1, 0});
            materials.push_back({0, 0, 1, 1});
            materials.push_back({2, 2, 2, 2});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 1);
            REQUIRE(mi.faces.size() == 4);
            REQUIRE(mi.vertices.size() == 4);
            validate(mi);
        }
        SECTION("Case 6: with duplicate materials")
        {
            materials.push_back({1, 1, 0, 0}); // material 4
            materials.push_back({0, 1, 1, 0}); // material 5
            materials.push_back({0, 1, 1, 0}); // material 6
            materials.push_back({1, 1, 0, 0}); // material 7
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 2);
            REQUIRE(mi.faces.size() == 7);
            REQUIRE(mi.vertices.size() == 5);
            validate(mi);
            REQUIRE(mi.unique_materials.size() == 6); // 4 bd materials + 2 input materials.
            REQUIRE(mi.unique_material_indices[4] == mi.unique_material_indices[7]);
            REQUIRE(mi.unique_material_indices[5] == mi.unique_material_indices[6]);
        }
    }
}
} // namespace

TEST_CASE("Material interface 2D", "[material_interface][2D]")
{
    using namespace simplicial_arrangement;
    SECTION("Int") { test_2D<Int>(); }
    SECTION("double") { test_2D<double>(); }
}

TEST_CASE("Material interface 3D", "[material_interface][3D]")
{
    using namespace simplicial_arrangement;
    SECTION("Int") { test_3D<Int>(); }
    SECTION("double") { test_3D<double>(); }
}
