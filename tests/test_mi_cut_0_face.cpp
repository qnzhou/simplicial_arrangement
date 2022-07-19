#include <simplicial_arrangement/material_interface.h>
#include <spdlog/spdlog.h>
#include <catch2/catch.hpp>

#include <implicit_predicates/implicit_predicates.h>

#include <initialize_simplicial_mi_complex.h>
#include <mi_cut_0_face.h>

namespace {

template <typename Scalar>
void test_2D()
{
    using namespace simplicial_arrangement;
    std::vector<Material<Scalar, 2>> materials;
    materials.push_back({0, 0, 0}); // material 3
    materials.push_back({1, 0, 0}); // material 4
    materials.push_back({0, 1, 0}); // material 5
    materials.push_back({0, 0, 1}); // material 6
    materials.push_back({1, 1, 1}); // material 7

    MaterialRepo<Scalar, 2> repo(materials);
    auto mi_complex = initialize_simplicial_mi_complex<2>();

    SECTION("Corners")
    {
        REQUIRE(mi_cut_0_face(repo, mi_complex, 0, 3) == 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, 0, 4) < 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, 0, 5) == 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, 0, 6) == 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, 0, 7) < 0);
    }

    SECTION("Edge mid point")
    {
        size_t vid = mi_complex.vertices.size();
        mi_complex.vertices.push_back({0, 5, 6});
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 3) > 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 4) > 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 5) == 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 6) == 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 7) < 0);
    }

    SECTION("Centroid")
    {
        size_t vid = mi_complex.vertices.size();
        mi_complex.vertices.push_back({4, 5, 6});
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 3) > 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 4) == 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 5) == 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 6) == 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 7) < 0);
    }
}

template <typename Scalar>
void test_3D()
{
    using namespace simplicial_arrangement;
    std::vector<Material<Scalar, 3>> materials;
    materials.push_back({0, 0, 0, 0}); // material 4
    materials.push_back({1, 0, 0, 0}); // material 5
    materials.push_back({0, 1, 0, 0}); // material 6
    materials.push_back({0, 0, 1, 0}); // material 7
    materials.push_back({0, 0, 0, 1}); // material 8
    materials.push_back({1, 1, 1, 1}); // material 9

    MaterialRepo<Scalar, 3> repo(materials);
    auto mi_complex = initialize_simplicial_mi_complex<3>();

    SECTION("Corners")
    {
        REQUIRE(mi_cut_0_face(repo, mi_complex, 0, 4) == 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, 0, 5) < 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, 0, 6) == 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, 0, 7) == 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, 0, 8) == 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, 0, 9) < 0);
    }

    SECTION("Edge mid point")
    {
        size_t vid = mi_complex.vertices.size();
        mi_complex.vertices.push_back({0, 1, 7, 8});
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 4) > 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 5) > 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 6) > 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 7) == 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 8) == 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 9) < 0);
    }

    SECTION("Face mid point")
    {
        size_t vid = mi_complex.vertices.size();
        mi_complex.vertices.push_back({0, 6, 7, 8});
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 4) > 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 5) > 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 6) == 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 7) == 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 8) == 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 9) < 0);
    }

    SECTION("Centroid")
    {
        size_t vid = mi_complex.vertices.size();
        mi_complex.vertices.push_back({5, 6, 7, 8});
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 4) > 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 5) == 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 6) == 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 7) == 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 8) == 0);
        REQUIRE(mi_cut_0_face(repo, mi_complex, vid, 9) < 0);
    }
}


} // namespace

TEST_CASE("mi_cut_0_face", "[material_interface]")
{
    using namespace simplicial_arrangement;
    SECTION("2D")
    {
#ifndef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
        SECTION("Int") { test_2D<Int>(); }
#endif
        SECTION("double") { test_2D<double>(); }
    }

    SECTION("3D")
    {
#ifndef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
        SECTION("Int") { test_3D<Int>(); }
#endif
        SECTION("double") { test_3D<double>(); }
    }
}
