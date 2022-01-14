#include <simplicial_arrangement/simplicial_arrangement.h>
#include <simplicial_arrangement/material_interface.h>
#include <simplicial_arrangement/lookup_table.h>

#include <catch2/catch.hpp>
#include <random>
#include <vector>

TEST_CASE("benchmark", "[arrangement][.benchmark]")
{
    using namespace simplicial_arrangement;
    load_lookup_table(BOTH);

    SECTION("2D")
    {
        constexpr size_t N = 1e2;
        std::vector<Plane<Int, 2>> int_data;
        std::vector<Plane<double, 2>> double_data;
        int_data.reserve(N);
        double_data.reserve(N);

        std::mt19937 gen(7);
        std::uniform_int_distribution<> distrib(-(1 << 12), 1 << 12);
        std::uniform_int_distribution<> rand_index(0, N-1);

        for (size_t i = 0; i < N; i++) {
            Int v0, v1, v2;
            v0 = distrib(gen);
            v1 = distrib(gen);
            v2 = distrib(gen);
            int_data.push_back({v0, v1, v2});
            double_data.push_back(
                {static_cast<double>(v0), static_cast<double>(v1), static_cast<double>(v2)});
        }

        BENCHMARK_ADVANCED("2D arrangement (int, 3 planes)")(Catch::Benchmark::Chronometer meter)
        {
            std::vector<Plane<Int, 2>> planes;
            planes.reserve(3);
            planes.push_back(int_data[rand_index(gen)]);
            planes.push_back(int_data[rand_index(gen)]);
            planes.push_back(int_data[rand_index(gen)]);

            meter.measure([&]() { return compute_arrangement(planes); });
        };
        BENCHMARK_ADVANCED("2D arrangement (double, 3 planes)")(Catch::Benchmark::Chronometer meter)
        {
            std::vector<Plane<double, 2>> planes;
            planes.reserve(3);
            planes.push_back(double_data[rand_index(gen)]);
            planes.push_back(double_data[rand_index(gen)]);
            planes.push_back(double_data[rand_index(gen)]);

            meter.measure([&]() { return compute_arrangement(planes); });
        };

        size_t num_vertices = 0;
        size_t num_faces = 0;
        size_t num_cells = 0;
        BENCHMARK("2D arrangement (int, 100 planes)")
        {
            auto r = compute_arrangement(int_data);
            num_vertices = r.vertices.size();
            num_faces = r.faces.size();
            num_cells = r.cells.size();
            return r;
        };
        BENCHMARK("2D arrangement (double, 100 planes)")
        {
            auto r = compute_arrangement(double_data);
            REQUIRE(r.vertices.size() == num_vertices);
            REQUIRE(r.faces.size() == num_faces);
            REQUIRE(r.cells.size() == num_cells);
            return r;
        };
    }

    SECTION("2D material interface")
    {
        constexpr size_t N = 1e2;
        std::vector<Material<Int, 2>> int_data;
        std::vector<Material<double, 2>> double_data;
        int_data.reserve(N);
        double_data.reserve(N);

        std::mt19937 gen(7);
        std::uniform_int_distribution<> distrib(-(1 << 12), 1 << 12);
        std::uniform_int_distribution<> rand_index(0, N-1);

        for (size_t i = 0; i < N; i++) {
            Int v0, v1, v2;
            v0 = distrib(gen);
            v1 = distrib(gen);
            v2 = distrib(gen);
            int_data.push_back({v0, v1, v2});
            double_data.push_back(
                {static_cast<double>(v0), static_cast<double>(v1), static_cast<double>(v2)});
        }

        BENCHMARK_ADVANCED("2D material interface (int, 3 planes)")(Catch::Benchmark::Chronometer meter)
        {
            std::vector<Material<Int, 2>> materials;
            materials.reserve(3);
            materials.push_back(int_data[rand_index(gen)]);
            materials.push_back(int_data[rand_index(gen)]);
            materials.push_back(int_data[rand_index(gen)]);

            meter.measure([&]() { return compute_material_interface(materials); });
        };
        BENCHMARK_ADVANCED("2D material interafce (double, 3 planes)")(Catch::Benchmark::Chronometer meter)
        {
            std::vector<Material<double, 2>> materials;
            materials.reserve(3);
            materials.push_back(double_data[rand_index(gen)]);
            materials.push_back(double_data[rand_index(gen)]);
            materials.push_back(double_data[rand_index(gen)]);

            meter.measure([&]() { return compute_material_interface(materials); });
        };

        size_t num_vertices = 0;
        size_t num_faces = 0;
        size_t num_cells = 0;
        BENCHMARK("2D material interface (int, 100 materials)")
        {
            auto r = compute_material_interface(int_data);
            num_vertices = r.vertices.size();
            num_faces = r.faces.size();
            num_cells = r.cells.size();
            return r;
        };
        BENCHMARK("2D material interface (double, 100 materials)")
        {
            auto r = compute_material_interface(double_data);
            REQUIRE(r.vertices.size() == num_vertices);
            REQUIRE(r.faces.size() == num_faces);
            REQUIRE(r.cells.size() == num_cells);
            return r;
        };
    }

    SECTION("3D")
    {
        constexpr size_t N = 1e2;
        std::vector<Plane<Int, 3>> int_data;
        std::vector<Plane<double, 3>> double_data;
        int_data.reserve(N);
        double_data.reserve(N);

        std::mt19937 gen(7);
        std::uniform_int_distribution<> distrib(-(1 << 12), 1 << 12);
        std::uniform_int_distribution<> rand_index(0, N-1);

        for (size_t i = 0; i < N; i++) {
            Int v0, v1, v2, v3;
            v0 = distrib(gen);
            v1 = distrib(gen);
            v2 = distrib(gen);
            v3 = distrib(gen);
            int_data.push_back({v0, v1, v2, v3});
            double_data.push_back({static_cast<double>(v0),
                static_cast<double>(v1),
                static_cast<double>(v2),
                static_cast<double>(v3)});
        }

        BENCHMARK_ADVANCED("3D arrangement (int, 1 planes, w/o lookup)")(Catch::Benchmark::Chronometer meter)
        {
            disable_lookup_table();
            std::vector<Plane<Int, 3>> planes;
            planes.reserve(1);
            planes.push_back(int_data[rand_index(gen)]);

            meter.measure([&]() { return compute_arrangement(planes); });
        };
        BENCHMARK_ADVANCED("3D arrangement (int, 1 planes), with lookup")(Catch::Benchmark::Chronometer meter)
        {
            enable_lookup_table();
            std::vector<Plane<Int, 3>> planes;
            planes.reserve(1);
            planes.push_back(int_data[rand_index(gen)]);

            meter.measure([&]() { return compute_arrangement(planes); });
        };

        BENCHMARK_ADVANCED("3D arrangement (double, 1 planes, w/o lookup)")(Catch::Benchmark::Chronometer meter)
        {
            disable_lookup_table();
            std::vector<Plane<double, 3>> planes;
            planes.reserve(1);
            planes.push_back(double_data[rand_index(gen)]);

            meter.measure([&]() { return compute_arrangement(planes); });
        };
        BENCHMARK_ADVANCED("3D arrangement (double, 1 planes, with lookup)")
        (Catch::Benchmark::Chronometer meter)
        {
            enable_lookup_table();
            std::vector<Plane<double, 3>> planes;
            planes.reserve(1);
            planes.push_back(double_data[rand_index(gen)]);

            meter.measure([&]() { return compute_arrangement(planes); });
        };

        BENCHMARK_ADVANCED("3D arrangement (int, 2 planes, w/o lookup)")(Catch::Benchmark::Chronometer meter)
        {
            disable_lookup_table();
            std::vector<Plane<Int, 3>> planes;
            planes.reserve(2);
            planes.push_back(int_data[rand_index(gen)]);
            planes.push_back(int_data[rand_index(gen)]);

            meter.measure([&]() { return compute_arrangement(planes); });
        };
        BENCHMARK_ADVANCED("3D arrangement (int, 2 planes, with lookup)")
        (Catch::Benchmark::Chronometer meter)
        {
            enable_lookup_table();
            std::vector<Plane<Int, 3>> planes;
            planes.reserve(2);
            planes.push_back(int_data[rand_index(gen)]);
            planes.push_back(int_data[rand_index(gen)]);

            meter.measure([&]() { return compute_arrangement(planes); });
        };

        BENCHMARK_ADVANCED("3D arrangement (double, 2 planes, w/o lookup)")(Catch::Benchmark::Chronometer meter)
        {
            disable_lookup_table();
            std::vector<Plane<double, 3>> planes;
            planes.reserve(2);
            planes.push_back(double_data[rand_index(gen)]);
            planes.push_back(double_data[rand_index(gen)]);

            meter.measure([&]() { return compute_arrangement(planes); });
        };
        BENCHMARK_ADVANCED("3D arrangement (double, 2 planes, with lookup)")
        (Catch::Benchmark::Chronometer meter)
        {
            enable_lookup_table();
            std::vector<Plane<double, 3>> planes;
            planes.reserve(2);
            planes.push_back(double_data[rand_index(gen)]);
            planes.push_back(double_data[rand_index(gen)]);

            meter.measure([&]() { return compute_arrangement(planes); });
        };

        BENCHMARK_ADVANCED("3D arrangement (int, 3 planes)")(Catch::Benchmark::Chronometer meter)
        {
            std::vector<Plane<Int, 3>> planes;
            planes.reserve(3);
            planes.push_back(int_data[rand_index(gen)]);
            planes.push_back(int_data[rand_index(gen)]);
            planes.push_back(int_data[rand_index(gen)]);

            meter.measure([&]() { return compute_arrangement(planes); });
        };
        BENCHMARK_ADVANCED("3D arrangement (double, 3 planes)")(Catch::Benchmark::Chronometer meter)
        {
            std::vector<Plane<double, 3>> planes;
            planes.reserve(3);
            planes.push_back(double_data[rand_index(gen)]);
            planes.push_back(double_data[rand_index(gen)]);
            planes.push_back(double_data[rand_index(gen)]);

            meter.measure([&]() { return compute_arrangement(planes); });
        };

        size_t num_vertices = 0;
        size_t num_faces = 0;
        size_t num_cells = 0;
        BENCHMARK("3D arrangement (int, 100 planes)")
        {
            auto r = compute_arrangement(int_data);
            num_vertices = r.vertices.size();
            num_faces = r.faces.size();
            num_cells = r.cells.size();
            return r;
        };
        BENCHMARK("3D arrangement (double, 100 planes)")
        {
            auto r = compute_arrangement(double_data);
            REQUIRE(r.vertices.size() == num_vertices);
            REQUIRE(r.faces.size() == num_faces);
            REQUIRE(r.cells.size() == num_cells);
            return r;
        };
    }

    SECTION("3D material interface")
    {
        constexpr size_t N = 1e2;
        std::vector<Material<Int, 3>> int_data;
        std::vector<Material<double, 3>> double_data;
        int_data.reserve(N);
        double_data.reserve(N);

        std::mt19937 gen(7);
        std::uniform_int_distribution<> distrib(-(1 << 12), 1 << 12);
        std::uniform_int_distribution<> rand_index(0, N-1);

        for (size_t i = 0; i < N; i++) {
            Int v0, v1, v2, v3;
            v0 = distrib(gen);
            v1 = distrib(gen);
            v2 = distrib(gen);
            v3 = distrib(gen);
            int_data.push_back({v0, v1, v2, v3});
            double_data.push_back({static_cast<double>(v0),
                static_cast<double>(v1),
                static_cast<double>(v2),
                static_cast<double>(v3)});
        }

        BENCHMARK_ADVANCED("3D material interface(int, 1 materials)")(Catch::Benchmark::Chronometer meter)
        {
            std::vector<Material<Int, 3>> materials;
            materials.reserve(1);
            materials.push_back(int_data[rand_index(gen)]);

            meter.measure([&]() { return compute_material_interface(materials); });
        };
        BENCHMARK_ADVANCED("3D material interface (double, 1 materials)")(Catch::Benchmark::Chronometer meter)
        {
            std::vector<Material<double, 3>> materials;
            materials.reserve(1);
            materials.push_back(double_data[rand_index(gen)]);

            meter.measure([&]() { return compute_material_interface(materials); });
        };

        BENCHMARK_ADVANCED("3D material interface (int, 2 materials, w/o lookup)")(Catch::Benchmark::Chronometer meter)
        {
            disable_lookup_table();
            std::vector<Material<Int, 3>> materials;
            materials.reserve(2);
            materials.push_back(int_data[rand_index(gen)]);
            materials.push_back(int_data[rand_index(gen)]);

            meter.measure([&]() { return compute_material_interface(materials); });
        };
        BENCHMARK_ADVANCED("3D material interface (int, 2 materials, with lookup)")(Catch::Benchmark::Chronometer meter)
        {
            enable_lookup_table();
            std::vector<Material<Int, 3>> materials;
            materials.reserve(2);
            materials.push_back(int_data[rand_index(gen)]);
            materials.push_back(int_data[rand_index(gen)]);

            meter.measure([&]() { return compute_material_interface(materials); });
        };
        BENCHMARK_ADVANCED("3D material interface (double, 2 materials, w/o lookup)")(Catch::Benchmark::Chronometer meter)
        {
            disable_lookup_table();
            std::vector<Material<double, 3>> materials;
            materials.reserve(2);
            materials.push_back(double_data[rand_index(gen)]);
            materials.push_back(double_data[rand_index(gen)]);

            meter.measure([&]() { return compute_material_interface(materials); });
        };
        BENCHMARK_ADVANCED("3D material interface (double, 2 materials, with lookup)")(Catch::Benchmark::Chronometer meter)
        {
            enable_lookup_table();
            std::vector<Material<double, 3>> materials;
            materials.reserve(2);
            materials.push_back(double_data[rand_index(gen)]);
            materials.push_back(double_data[rand_index(gen)]);

            meter.measure([&]() { return compute_material_interface(materials); });
        };

        BENCHMARK_ADVANCED("3D material interface (int, 3 materials, w/o lookup)")(Catch::Benchmark::Chronometer meter)
        {
            disable_lookup_table();
            std::vector<Material<Int, 3>> materials;
            materials.reserve(3);
            materials.push_back(int_data[rand_index(gen)]);
            materials.push_back(int_data[rand_index(gen)]);
            materials.push_back(int_data[rand_index(gen)]);

            meter.measure([&]() { return compute_material_interface(materials); });
        };
        BENCHMARK_ADVANCED("3D material interface (int, 3 materials, with lookup)")(Catch::Benchmark::Chronometer meter)
        {
            enable_lookup_table();
            std::vector<Material<Int, 3>> materials;
            materials.reserve(3);
            materials.push_back(int_data[rand_index(gen)]);
            materials.push_back(int_data[rand_index(gen)]);
            materials.push_back(int_data[rand_index(gen)]);

            meter.measure([&]() { return compute_material_interface(materials); });
        };
        BENCHMARK_ADVANCED("3D material interface (double, 3 materials, w/o lookup)")(Catch::Benchmark::Chronometer meter)
        {
            disable_lookup_table();
            std::vector<Material<double, 3>> materials;
            materials.reserve(3);
            materials.push_back(double_data[rand_index(gen)]);
            materials.push_back(double_data[rand_index(gen)]);
            materials.push_back(double_data[rand_index(gen)]);

            meter.measure([&]() { return compute_material_interface(materials); });
        };
        BENCHMARK_ADVANCED("3D material interface (double, 3 materials, with lookup)")(Catch::Benchmark::Chronometer meter)
        {
            enable_lookup_table();
            std::vector<Material<double, 3>> materials;
            materials.reserve(3);
            materials.push_back(double_data[rand_index(gen)]);
            materials.push_back(double_data[rand_index(gen)]);
            materials.push_back(double_data[rand_index(gen)]);

            meter.measure([&]() { return compute_material_interface(materials); });
        };

        size_t num_vertices = 0;
        size_t num_faces = 0;
        size_t num_cells = 0;
        BENCHMARK("3D material interface (int, 100 planes)")
        {
            auto r = compute_material_interface(int_data);
            num_vertices = r.vertices.size();
            num_faces = r.faces.size();
            num_cells = r.cells.size();
            return r;
        };
        BENCHMARK("3D material interface (double, 100 planes)")
        {
            auto r = compute_material_interface(double_data);
            REQUIRE(r.vertices.size() == num_vertices);
            REQUIRE(r.faces.size() == num_faces);
            REQUIRE(r.cells.size() == num_cells);
            return r;
        };
    }
}
