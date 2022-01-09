#include <simplicial_arrangement/material_interface.h>

#include <implicit_predicates/implicit_predicates.h>

#include <absl/container/flat_hash_set.h>
#include <spdlog/spdlog.h>
#include <catch2/catch.hpp>

namespace {

template <int DIM>
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

    if constexpr (DIM == 2) {
        // Check face orientation is consistent within each cell.
        for (const auto& c : mi.cells) {
            const size_t num_cell_faces = c.faces.size();
            for (size_t i = 0; i < num_cell_faces; i++) {
                const auto& f0 = mi.faces[c.faces[i]];
                const auto& f1 = mi.faces[c.faces[(i + 1) % num_cell_faces]];
                if (f0.negative_material_label == c.material_label) {
                    REQUIRE((f0.vertices[1] == f1.vertices[0] || f0.vertices[1] == f1.vertices[1]));
                } else {
                    REQUIRE((f0.vertices[0] == f1.vertices[0] || f0.vertices[0] == f1.vertices[1]));
                }
            }
        }
    } else {
        // Check each cell is topologically a sphere using Euler
        // characteristics.
        absl::flat_hash_set<size_t> vertex_set;
        absl::flat_hash_set<std::array<size_t, 2>> edge_set;
        vertex_set.reserve(mi.vertices.size());
        edge_set.reserve(mi.vertices.size() * 2);
        for (const auto& c : mi.cells) {
            vertex_set.clear();
            edge_set.clear();
            for (const auto fid : c.faces) {
                const auto& f = mi.faces[fid];
                vertex_set.insert(f.vertices.begin(), f.vertices.end());
                const size_t num_face_vertices = f.vertices.size();
                for (size_t i = 0; i < num_face_vertices; i++) {
                    const size_t curr_vid = f.vertices[i];
                    const size_t next_vid = f.vertices[(i + 1) % num_face_vertices];
                    std::array<size_t, 2> e{curr_vid, next_vid};
                    if (e[0] > e[1]) std::swap(e[0], e[1]);
                    edge_set.insert(std::move(e));
                }
            }

            int euler = int(vertex_set.size()) - int(edge_set.size()) + int(c.faces.size());
            REQUIRE(euler == 2);
        }
    }
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
    SECTION("3 implicits")
    {
        SECTION("Case 1")
        {
            materials.push_back({1, 0, 0, 0});
            materials.push_back({0, 1, 0, 0});
            materials.push_back({0, 0, 1, 0});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 3);
            REQUIRE(mi.faces.size() == 12);
            REQUIRE(mi.vertices.size() == 8);
            validate(mi);
        }
        SECTION("Case 2")
        {
            materials.push_back({1, 0, 0, 0});
            materials.push_back({0, 1, 0, 0});
            materials.push_back({2, 0, 0, 0});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 2);
            REQUIRE(mi.faces.size() == 7);
            REQUIRE(mi.vertices.size() == 5);
            validate(mi);
        }
        SECTION("Case 3")
        {
            materials.push_back({1, 0, 0, 0});
            materials.push_back({0, 1, 0, 0});
            materials.push_back({0, 0, 1, 1});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 3);
            REQUIRE(mi.faces.size() == 13);
            REQUIRE(mi.vertices.size() == 11);
            validate(mi);
        }
        SECTION("Case 4")
        {
            materials.push_back({1, 0, 0, 2});
            materials.push_back({0, 2, 2, 1});
            materials.push_back({2, 1, 1, 0});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 3);
            REQUIRE(mi.faces.size() == 13);
            REQUIRE(mi.vertices.size() == 11);
            validate(mi);
        }
    }
    SECTION("4 implicits")
    {
        SECTION("Case 1")
        {
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
        SECTION("Case 2")
        {
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
        SECTION("Case 3")
        {
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
        SECTION("Case 4")
        {
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
        SECTION("Case 5")
        {
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
