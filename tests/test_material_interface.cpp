#include <simplicial_arrangement/lookup_table.h>
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
        SECTION("Case 5")
        {
            materials.push_back({1, 1, 1, 2});
            materials.push_back({2, 2, 2, 1});
            materials.push_back({0, 0, 0, 0});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 2);
            REQUIRE(mi.faces.size() == 8);
            REQUIRE(mi.vertices.size() == 7);
            validate(mi);
        }
        SECTION("Case 6")
        {
            materials.push_back({1, 1, 1, 1});
            materials.push_back({0, 0, 2, 0});
            materials.push_back({2, 0, 0, 0});
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 3);
            REQUIRE(mi.faces.size() == 12);
            REQUIRE(mi.vertices.size() == 9);
            validate(mi);

            // Check vertex shared by more than 4 materials uses the materials
            // with the smallest index.
            const auto& vertices = mi.vertices;
            const auto& faces = mi.faces;
            std::vector<size_t> vertex_valence(vertices.size(), 0);
            for (const auto& f : faces) {
                for (const auto& vid : f.vertices) {
                    vertex_valence[vid]++;
                }
            }
            const auto itr = std::max_element(vertex_valence.begin(), vertex_valence.end());
            REQUIRE(*itr == 8);
            auto v = vertices[itr - vertex_valence.begin()];
            std::sort(v.begin(), v.end());
            REQUIRE(v[0] == 1);
            REQUIRE(v[1] == 3);
            REQUIRE(v[2] == 4);
            REQUIRE(v[3] == 5);
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
        SECTION("Case 7")
        {
            materials.push_back({0, 0, 2, 0}); // material 4
            materials.push_back({2, 0, 0, 0}); // material 5
            materials.push_back({4, -1, -2, -1}); // material 6
            materials.push_back({-2, -1, 4, -1}); // material 7
            auto mi = compute_material_interface(materials);
            REQUIRE(mi.cells.size() == 4);
            REQUIRE(mi.faces.size() == 15);
            REQUIRE(mi.vertices.size() == 9);
            validate(mi);
        }
    }
}
} // namespace

TEST_CASE("Material interface 2D", "[material_interface][2D]")
{
    using namespace simplicial_arrangement;
#ifndef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
    SECTION("Int") { test_2D<Int>(); }
#endif
    SECTION("double") { test_2D<double>(); }
}

TEST_CASE("Material interface 3D", "[material_interface][3D]")
{
    using namespace simplicial_arrangement;
    REQUIRE(load_lookup_table(MATERIAL_INTERFACE));

#ifndef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
    SECTION("Int")
    {
        SECTION("Without lookup")
        {
            disable_lookup_table();
            test_3D<Int>();
        }
        SECTION("With lookup")
        {
            enable_lookup_table();
            test_3D<Int>();
        }
        disable_lookup_table();
    }
#endif
    SECTION("double")
    {
        SECTION("Without lookup")
        {
            disable_lookup_table();
            test_3D<double>();
        }
        SECTION("With lookup")
        {
            enable_lookup_table();
            test_3D<double>();
        }
        disable_lookup_table();
    }
}

TEST_CASE("Material interface degeneracy", "[mi][degeneracy][3D]")
{
    using namespace simplicial_arrangement;
    disable_lookup_table();

    using Scalar = double;
    std::vector<Material<Scalar, 3>> materials;

    auto assert_equivalent = [](const auto& mi1, const auto& mi2) {
        REQUIRE(mi1.vertices.size() == mi2.vertices.size());
        REQUIRE(mi1.faces.size() == mi2.faces.size());
        REQUIRE(mi1.cells.size() == mi2.cells.size());
    };

    SECTION("5 materials in different order")
    {
        SECTION("Point intersection")
        {
            materials.push_back(
                {-1.155765120458494, -3.192762224028497, 0.5553418458787984, 3.793185498608192});
            materials.push_back(
                {-2.392200834775261, -1.528746399924049, 5.208044446314545, -1.287097211615235});
            materials.push_back(
                {1.271220156626018, -1.439963072432723, 0.6169791304579846, -0.4482362146512795});
            materials.push_back(
                {6.390380078685507, 5.280405125348178, 4.520264151726717, -16.1910493557604});
            materials.push_back(
                {-6.781799097287127, 5.240177899219642, -8.129031919657647, 9.670653117725132});

            const auto r1 = compute_material_interface(materials);
            std::reverse(materials.begin(), materials.end());
            const auto r2 = compute_material_interface(materials);
            assert_equivalent(r1, r2);
        }
        SECTION("Segment intersection")
        {
            materials.push_back(
                {7.718609961072335, -9.204820073642072, -4.746189735932862, 6.232399848502599});
            materials.push_back(
                {-1.597335209387186, -0.4805352955931301, 5.753076219347818, -3.675205714367502});
            materials.push_back(
                {-8.454966357160966, 13.00990946531108, -0.6549198591392553, -3.900023249010856});
            materials.push_back(
                {-2.677583008334603, -0.3260814570562367, 8.684911939116279, -5.681247473725441});
            materials.push_back(
                {4.425259585475192, -4.509084694383795, -4.257609367657988, 4.34143447656659});

            const auto r1 = compute_material_interface(materials);
            std::reverse(materials.begin(), materials.end());
            const auto r2 = compute_material_interface(materials);
            assert_equivalent(r1, r2);
        }
        SECTION("Triangle intersection")
        {
            materials.push_back(
                {2.137454324991177, -2.137454315051372, -2.137454319214663, -2.137454327568319});
            materials.push_back(
                {-9.391832302908245, 9.391832300523651, 9.391832312300817, 9.391832301038432});
            materials.push_back(
                {1.681410363741283, -1.681410373407355, -1.681410377891855, -1.681410372213574});
            materials.push_back(
                {4.867248900161551, -4.867248902720028, -4.867248888583141, -4.867248900564688});
            materials.push_back(
                {7.602325749023126, -7.60232575169791, -7.602325755232684, -7.602325750807439});

            const auto r1 = compute_material_interface(materials);
            std::reverse(materials.begin(), materials.end());
            const auto r2 = compute_material_interface(materials);
            assert_equivalent(r1, r2);
        }
        SECTION("Quad intersection")
        {
            materials.push_back(
                {-7.395795332281445, -7.395795334049828, 7.39579533993123, 7.3957953347206});
            materials.push_back(
                {6.411832922871854, 6.411832913469024, -6.411832914570133, -6.411832919619761});
            materials.push_back(
                {1.036045238227753, 1.036045244464811, -1.036045238074583, -1.036045253607837});
            materials.push_back(
                {-5.986464634521012, -5.986464625816754, 5.986464619923138, 5.986464624311083});
            materials.push_back(
                {-3.737885309931638, -3.737885301409309, 3.737885298445151, 3.737885295185574});

            const auto r1 = compute_material_interface(materials);
            std::reverse(materials.begin(), materials.end());
            const auto r2 = compute_material_interface(materials);
            assert_equivalent(r1, r2);
        }
    }
}
