#include <simplicial_arrangement/lookup_table.h>
#include <simplicial_arrangement/material_interface.h>
#include <simplicial_arrangement/simplicial_arrangement.h>

#include <implicit_predicates/implicit_predicates.h>

#include <catch2/catch.hpp>
#include <iomanip>
#include <iostream>
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
        size_t type2_failure = 0; // Failure due to invalid state.
        for (size_t i = 0; i < N; i++) {
            planes.clear();
            for (size_t j = 0; j < M; j++) {
                Scalar v0 = distrib(gen);
                Scalar v1 = distrib(gen);
                Scalar v2 = distrib(gen);
                Scalar v3 = -v0 - v1 - v2;
                planes.push_back({v0, v1, v2, v3});
            }
            try {
                const auto r1 = compute_arrangement(planes);

                std::reverse(planes.begin(), planes.end());
                const auto r2 = compute_arrangement(planes);

                // A necessary condition for success.
                if (r1.cells.size() == r2.cells.size() && r1.faces.size() == r2.faces.size() &&
                    r1.vertices.size() == r2.vertices.size()) {
                    success++;
                }
            } catch (const std::runtime_error&) {
                type2_failure++;
            }
        }

        INFO("Type 1 failure: " << N - success - type2_failure);
        INFO("Type 2 failure: " << type2_failure);
        REQUIRE(success == N);
    }

    SECTION("Simplicial arrangement - segment intersection")
    {
        size_t success = 0;
        size_t type2_failure = 0; // Failure due to invalid state.
        for (size_t i = 0; i < N; i++) {
            planes.clear();
            for (size_t j = 0; j < M; j++) {
                Scalar a = distrib(gen);
                Scalar b = distrib(gen);
                planes.push_back({a + 2 * b, -2 * a - 3 * b, a, b});
            }
            try {
                const auto r1 = compute_arrangement(planes);

                std::reverse(planes.begin(), planes.end());
                const auto r2 = compute_arrangement(planes);

                // A necessary condition for success.
                if (r1.cells.size() == r2.cells.size() && r1.faces.size() == r2.faces.size() &&
                    r1.vertices.size() == r2.vertices.size()) {
                    success++;
                }
            } catch (const std::runtime_error&) {
                type2_failure++;
            }
        }

        INFO("Type 1 failure: " << N - success - type2_failure);
        INFO("Type 2 failure: " << type2_failure);
        REQUIRE(success == N);
    }

    SECTION("Simplicial arrangement - triangle intersection")
    {
        size_t success = 0;
        size_t type2_failure = 0; // Failure due to invalid state.
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
            try {
                const auto r1 = compute_arrangement(planes);

                std::reverse(planes.begin(), planes.end());
                const auto r2 = compute_arrangement(planes);

                // A necessary condition for success.
                if (r1.cells.size() == r2.cells.size() && r1.faces.size() == r2.faces.size() &&
                    r1.vertices.size() == r2.vertices.size()) {
                    success++;
                }
            } catch (const std::runtime_error&) {
                type2_failure++;
            }
        }

        INFO("Type 1 failure: " << N - success - type2_failure);
        INFO("Type 2 failure: " << type2_failure);
        REQUIRE(success == N);
    }

    SECTION("Simplicial arrangement - quad intersection")
    {
        size_t success = 0;
        size_t type2_failure = 0; // Failure due to invalid state.
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
            try {
                const auto r1 = compute_arrangement(planes);

                std::reverse(planes.begin(), planes.end());
                const auto r2 = compute_arrangement(planes);

                // A necessary condition for success.
                if (r1.cells.size() == r2.cells.size() && r1.faces.size() == r2.faces.size() &&
                    r1.vertices.size() == r2.vertices.size()) {
                    success++;
                }
            } catch (const std::runtime_error&) {
                type2_failure++;
            }
        }

        INFO("Type 1 failure: " << N - success - type2_failure);
        INFO("Type 2 failure: " << type2_failure);
        REQUIRE(success == N);
    }

    SECTION("Simplicial arrangement - flat")
    {
        size_t success = 0;
        size_t type2_failure = 0; // Failure due to invalid state.
        for (size_t i = 0; i < N; i++) {
            planes.clear();
            for (size_t j = 0; j < M; j++) {
                // Tiny values at the tet vertices.
                Scalar t0 = distrib(gen) / (1 << 30);
                Scalar t1 = distrib(gen) / (1 << 30);
                Scalar t2 = distrib(gen) / (1 << 30);
                Scalar t3 = distrib(gen) / (1 << 30);
                planes.push_back({t0, t1, t2, t3});
            }
            try {
                const auto r1 = compute_arrangement(planes);

                std::reverse(planes.begin(), planes.end());
                const auto r2 = compute_arrangement(planes);

                // A necessary condition for success.
                if (r1.cells.size() == r2.cells.size() && r1.faces.size() == r2.faces.size() &&
                    r1.vertices.size() == r2.vertices.size()) {
                    success++;
                }
            } catch (const std::runtime_error&) {
                type2_failure++;
            }
        }

        INFO("Type 1 failure: " << N - success - type2_failure);
        INFO("Type 2 failure: " << type2_failure);
        REQUIRE(success == N);
    }

    std::vector<Material<Scalar, 3>> materials;
    materials.reserve(M + 1);
    SECTION("Material interface - point intersection")
    {
        size_t success = 0;
        size_t type2_failure = 0; // Failure due to invalid state.
        for (size_t i = 0; i < N; i++) {
            materials.clear();
            for (size_t j = 0; j <= M; j++) {
                Scalar v0 = distrib(gen);
                Scalar v1 = distrib(gen);
                Scalar v2 = distrib(gen);
                Scalar v3 = -v0 - v1 - v2;
                materials.push_back({v0, v1, v2, v3});
            }
            try {
                const auto r1 = compute_material_interface(materials);

                std::reverse(materials.begin(), materials.end());
                const auto r2 = compute_material_interface(materials);

                // A necessary condition for success.
                if (r1.cells.size() == r2.cells.size() && r1.faces.size() == r2.faces.size() &&
                    r1.vertices.size() == r2.vertices.size()) {
                    success++;
                }
            } catch (const std::runtime_error&) {
                type2_failure++;
            }
        }

        INFO("Type 1 failure: " << N - success - type2_failure);
        INFO("Type 2 failure: " << type2_failure);
        REQUIRE(success == N);
    }

    SECTION("Material interface - segment intersection")
    {
        size_t success = 0;
        size_t type2_failure = 0; // Failure due to invalid state.
        for (size_t i = 0; i < N; i++) {
            materials.clear();
            for (size_t j = 0; j <= M; j++) {
                Scalar a = distrib(gen);
                Scalar b = distrib(gen);
                materials.push_back({a + 2 * b, -2 * a - 3 * b, a, b});
            }
            try {
                const auto r1 = compute_material_interface(materials);

                std::reverse(materials.begin(), materials.end());
                const auto r2 = compute_material_interface(materials);

                // A necessary condition for success.
                if (r1.cells.size() == r2.cells.size() && r1.faces.size() == r2.faces.size() &&
                    r1.vertices.size() == r2.vertices.size()) {
                    success++;
                }
            } catch (const std::runtime_error&) {
                type2_failure++;
            }
        }

        INFO("Type 1 failure: " << N - success - type2_failure);
        INFO("Type 2 failure: " << type2_failure);
        REQUIRE(success == N);
    }

    SECTION("Material interface - triangle intersection")
    {
        size_t success = 0;
        size_t type2_failure = 0; // Failure due to invalid state.
        for (size_t i = 0; i < N; i++) {
            materials.clear();
            for (size_t j = 0; j <= M; j++) {
                Scalar a = distrib(gen) * (j % 2 == 0 ? 1 : -1);
                // Small perturbations per vertex.
                Scalar t0 = distrib(gen) / (1 << 30);
                Scalar t1 = distrib(gen) / (1 << 30);
                Scalar t2 = distrib(gen) / (1 << 30);
                Scalar t3 = distrib(gen) / (1 << 30);
                materials.push_back({a + t0, -a + t1, -a + t2, -a + t3});
            }
            try {
                const auto r1 = compute_material_interface(materials);

                std::reverse(materials.begin(), materials.end());
                const auto r2 = compute_material_interface(materials);

                // A necessary condition for success.
                if (r1.cells.size() == r2.cells.size() && r1.faces.size() == r2.faces.size() &&
                    r1.vertices.size() == r2.vertices.size()) {
                    success++;
                }
            } catch (const std::runtime_error&) {
                type2_failure++;
            }
        }

        INFO("Type 1 failure: " << N - success - type2_failure);
        INFO("Type 2 failure: " << type2_failure);
        REQUIRE(success == N);
    }

    SECTION("Material interface - quad intersection")
    {
        size_t success = 0;
        size_t type2_failure = 0; // Failure due to invalid state.
        for (size_t i = 0; i < N; i++) {
            materials.clear();
            for (size_t j = 0; j <= M; j++) {
                Scalar a = distrib(gen) * (j % 2 == 0 ? 1 : -1);
                // Small perturbations per vertex.
                Scalar t0 = distrib(gen) / (1 << 30);
                Scalar t1 = distrib(gen) / (1 << 30);
                Scalar t2 = distrib(gen) / (1 << 30);
                Scalar t3 = distrib(gen) / (1 << 30);
                materials.push_back({a + t0, a + t1, -a + t2, -a + t3});
            }
            try {
                const auto r1 = compute_material_interface(materials);

                std::reverse(materials.begin(), materials.end());
                const auto r2 = compute_material_interface(materials);

                // A necessary condition for success.
                if (r1.cells.size() == r2.cells.size() && r1.faces.size() == r2.faces.size() &&
                    r1.vertices.size() == r2.vertices.size()) {
                    success++;
                }
            } catch (const std::runtime_error&) {
                type2_failure++;
            }
        }

        INFO("Type 1 failure: " << N - success - type2_failure);
        INFO("Type 2 failure: " << type2_failure);
        REQUIRE(success == N);
    }

    SECTION("Material interface - flat")
    {
        std::cout << std::numeric_limits<Scalar>::epsilon() << std::endl;
        size_t success = 0;
        size_t type2_failure = 0; // Failure due to invalid state.
        for (size_t i = 0; i < N; i++) {
            materials.clear();
            for (size_t j = 0; j <= M; j++) {
                // Tiny values at the tet vertices.
                Scalar t0 = distrib(gen) * std::numeric_limits<Scalar>::epsilon();
                Scalar t1 = distrib(gen) * std::numeric_limits<Scalar>::epsilon();
                Scalar t2 = distrib(gen) * std::numeric_limits<Scalar>::epsilon();
                Scalar t3 = distrib(gen) * std::numeric_limits<Scalar>::epsilon();
                materials.push_back({t0, t1, t2, t3});
            }
            try {
                const auto r1 = compute_material_interface(materials);

                std::reverse(materials.begin(), materials.end());
                const auto r2 = compute_material_interface(materials);

                // A necessary condition for success.
                if (r1.cells.size() == r2.cells.size() && r1.faces.size() == r2.faces.size() &&
                    r1.vertices.size() == r2.vertices.size()) {
                    std::cout << r1.cells.size() << " " << r1.faces.size() << " "
                              << r1.vertices.size() << std::endl;
                    success++;
                }
            } catch (const std::runtime_error&) {
                type2_failure++;
            }
        }

        INFO("Type 1 failure: " << N - success - type2_failure);
        INFO("Type 2 failure: " << type2_failure);
        REQUIRE(success == N);
    }

    const auto stage_stats = implicit_predicates::get_stage_stats();
    std::cout << "Semi-static filter: " << stage_stats[0] << std::endl;
    std::cout << "Interval arithmetic: " << stage_stats[1] << std::endl;
    std::cout << "Exact computation: " << stage_stats[2] << std::endl;
}
