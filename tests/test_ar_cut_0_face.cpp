#include <simplicial_arrangement/simplicial_arrangement.h>
#include <spdlog/spdlog.h>
#include <catch2/catch.hpp>

#include <implicit_predicates/implicit_predicates.h>

#include <ar_cut_0_face.h>
#include <initialize_simplicial_ar_complex.h>

namespace {

template <typename Scalar>
void test_2D()
{
    using namespace simplicial_arrangement;
    std::vector<Plane<Scalar, 2>> planes;
    planes.push_back({0, 0, 0}); // plane 3
    planes.push_back({1, 0, 0}); // plane 4
    planes.push_back({0, 1, 0}); // plane 5
    planes.push_back({0, 0, 1}); // plane 6
    planes.push_back({1, 1, 1}); // plane 7
    planes.push_back({-1, -1, 2}); // plane 8

    PlaneRepo<Scalar, 2> repo(planes);
    auto ar_complex = initialize_simplicial_ar_complex<2>();

    REQUIRE(ar_cut_0_face(repo, ar_complex, 0, 3) == 0);
    REQUIRE(ar_cut_0_face(repo, ar_complex, 0, 4) > 0);
    REQUIRE(ar_cut_0_face(repo, ar_complex, 0, 5) == 0);
    REQUIRE(ar_cut_0_face(repo, ar_complex, 0, 6) == 0);
    REQUIRE(ar_cut_0_face(repo, ar_complex, 0, 7) > 0);
    REQUIRE(ar_cut_0_face(repo, ar_complex, 0, 8) < 0);
}

template <typename Scalar>
void test_3D()
{
    using namespace simplicial_arrangement;
    std::vector<Plane<Scalar, 3>> planes;
    planes.push_back({0, 0, 0, 0}); // plane 4
    planes.push_back({1, 0, 0, 0}); // plane 5
    planes.push_back({0, 1, 0, 0}); // plane 6
    planes.push_back({0, 0, 1, 0}); // plane 7
    planes.push_back({0, 0, 0, 1}); // plane 8
    planes.push_back({1, 1, 1, 1}); // plane 9
    planes.push_back({-1, -1, 2, 2}); // plane 10

    PlaneRepo<Scalar, 3> repo(planes);
    auto ar_complex = initialize_simplicial_ar_complex<3>();

    REQUIRE(ar_cut_0_face(repo, ar_complex, 0, 4) == 0);
    REQUIRE(ar_cut_0_face(repo, ar_complex, 0, 5) > 0);
    REQUIRE(ar_cut_0_face(repo, ar_complex, 0, 6) == 0);
    REQUIRE(ar_cut_0_face(repo, ar_complex, 0, 7) == 0);
    REQUIRE(ar_cut_0_face(repo, ar_complex, 0, 8) == 0);
    REQUIRE(ar_cut_0_face(repo, ar_complex, 0, 9) > 0);
    REQUIRE(ar_cut_0_face(repo, ar_complex, 0, 10) < 0);
}

} // namespace

TEST_CASE("ar_cut_0_face", "[arrangement]")
{
    using namespace simplicial_arrangement;
    SECTION("2D")
    {
        SECTION("Int") { test_2D<Int>(); }
        SECTION("double") { test_2D<double>(); }
    }

    SECTION("3D")
    {
        SECTION("Int") { test_3D<Int>(); }
        SECTION("double") { test_3D<double>(); }
    }
}
