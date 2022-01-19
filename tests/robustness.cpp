#include <simplicial_arrangement/lookup_table.h>
#include <simplicial_arrangement/material_interface.h>
#include <simplicial_arrangement/simplicial_arrangement.h>

#include <catch2/catch.hpp>
#include <random>
#include <vector>

/**
 * The following test check for robustness of the arrangement and material
 * interface algorithm.  It should only pass when
 * `SIMPLICIAL_ARRANGEMENT_NON_ROBUST` flat is OFF.
 */
TEST_CASE("robustness", "[ar][mi][.robustness]")
{
    using namespace simplicial_arrangement;
    disable_lookup_table();

    using Scalar = double;
    constexpr size_t N = 1e4;
    constexpr size_t M = 4;
    std::vector<Plane<Scalar, 3>> planes;
    planes.reserve(M);

    std::mt19937 gen(1);
    std::uniform_real_distribution<> distrib(-10, 10);

    SECTION("Simplicial arrangement - point intersection")
    {
        size_t success = 0;
        for (size_t i = 0; i < N; i++) {
            planes.clear();
            for (size_t j = 0; j < M; j++) {
                Scalar v0 = distrib(gen);
                Scalar v1 = distrib(gen);
                Scalar v2 = distrib(gen);
                Scalar v3 = -v0 - v1 - v2;
                planes.push_back({v0, v1, v2, v3});
            }
            const auto r1 = compute_arrangement(planes);

            std::reverse(planes.begin(), planes.end());
            const auto r2 = compute_arrangement(planes);

            // A necessary condition for success.
            if (r1.cells.size() == r2.cells.size() && r1.faces.size() == r2.faces.size() &&
                r1.vertices.size() == r2.vertices.size()) {
                success++;
            }
        }

        REQUIRE(success == N);
    }

    SECTION("Simplicial arrangement - near point intersection")
    {
        size_t success = 0;
        for (size_t i = 0; i < N; i++) {
            planes.clear();
            for (size_t j = 0; j < M; j++) {
                Scalar v0 = distrib(gen);
                Scalar v1 = distrib(gen);
                Scalar v2 = distrib(gen);
                Scalar v3 = distrib(gen);

                // The following operations produce floating point error.
                Scalar s = (v0 + v1 + v2 + v3) / 4;
                v0 -= s;
                v1 -= s;
                v2 -= s;
                v3 -= s;

                planes.push_back({v0, v1, v2, v3});
            }
            const auto r1 = compute_arrangement(planes);

            std::reverse(planes.begin(), planes.end());
            const auto r2 = compute_arrangement(planes);

            // A necessary condition for success.
            if (r1.cells.size() == r2.cells.size() && r1.faces.size() == r2.faces.size() &&
                r1.vertices.size() == r2.vertices.size()) {
                success++;
            }
        }

        REQUIRE(success == N);
    }


    SECTION("Simplicial arrangement - line intersection")
    {
        // Sadly, the non-robust version enters an infinite loop for this test.
#ifndef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
        size_t success = 0;
        for (size_t i = 0; i < N; i++) {
            planes.clear();
            for (size_t j = 0; j < M; j++) {
                Scalar a = distrib(gen);
                Scalar b = distrib(gen);
                planes.push_back({a + 2 * b, -2 * a - 3 * b, a, b});
            }
            const auto r1 = compute_arrangement(planes);

            std::reverse(planes.begin(), planes.end());
            const auto r2 = compute_arrangement(planes);

            // A necessary condition for success.
            if (r1.cells.size() == r2.cells.size() && r1.faces.size() == r2.faces.size() &&
                r1.vertices.size() == r2.vertices.size()) {
                success++;
            }
        }

        REQUIRE(success == N);
#endif
    }

    SECTION("Simplicial arrangement - nearly coplanar quads")
    {
        // Sadly, the non-robust version enters an infinite loop for this test.
#ifndef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
        size_t success = 0;
        for (size_t i = 0; i < N; i++) {
            planes.clear();
            for (size_t j = 0; j < M; j++) {
                Scalar a = distrib(gen);
                // Small perturbations per vertex.
                Scalar t0 = distrib(gen) / (1 << 30);
                Scalar t1 = distrib(gen) / (1 << 30);
                Scalar t2 = distrib(gen) / (1 << 30);
                Scalar t3 = distrib(gen) / (1 << 30);
                planes.push_back({a + t0, a + t1, -a + t2, -a + t3});
            }
            const auto r1 = compute_arrangement(planes);

            std::reverse(planes.begin(), planes.end());
            const auto r2 = compute_arrangement(planes);

            // A necessary condition for success.
            if (r1.cells.size() == r2.cells.size() && r1.faces.size() == r2.faces.size() &&
                r1.vertices.size() == r2.vertices.size()) {
                success++;
            }
        }

        REQUIRE(success == N);
#endif
    }

    SECTION("Simplicial arrangement - nearly coplanar triangles")
    {
        // Sadly, the non-robust version enters an infinite loop for this test.
#ifndef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
        size_t success = 0;
        for (size_t i = 0; i < N; i++) {
            planes.clear();
            for (size_t j = 0; j < M; j++) {
                Scalar a = distrib(gen);
                // Small perturbations per vertex.
                Scalar t0 = distrib(gen) / (1 << 30);
                Scalar t1 = distrib(gen) / (1 << 30);
                Scalar t2 = distrib(gen) / (1 << 30);
                Scalar t3 = distrib(gen) / (1 << 30);
                planes.push_back({a + t0, -a + t1, -a + t2, -a + t3});
            }
            const auto r1 = compute_arrangement(planes);

            std::reverse(planes.begin(), planes.end());
            const auto r2 = compute_arrangement(planes);

            // A necessary condition for success.
            if (r1.cells.size() == r2.cells.size() && r1.faces.size() == r2.faces.size() &&
                r1.vertices.size() == r2.vertices.size()) {
                success++;
            }
        }

        REQUIRE(success == N);
#endif
    }

    std::vector<Material<Scalar, 3>> materials;
    materials.reserve(M + 1);
    SECTION("Material interface - point intersection")
    {
        size_t success = 0;
        for (size_t i = 0; i < N; i++) {
            materials.clear();
            for (size_t j = 0; j <= M; j++) {
                Scalar v0 = distrib(gen);
                Scalar v1 = distrib(gen);
                Scalar v2 = distrib(gen);
                Scalar v3 = -v0 - v1 - v2;
                materials.push_back({v0, v1, v2, v3});
            }
            const auto r1 = compute_material_interface(materials);

            std::reverse(materials.begin(), materials.end());
            const auto r2 = compute_material_interface(materials);

            // A necessary condition for success.
            if (r1.cells.size() == r2.cells.size() && r1.faces.size() == r2.faces.size() &&
                r1.vertices.size() == r2.vertices.size()) {
                success++;
            }
        }

        REQUIRE(success == N);
    }
}
