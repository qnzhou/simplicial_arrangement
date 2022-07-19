#include <simplicial_arrangement/lookup_table.h>
#include <simplicial_arrangement/material_interface.h>
#include <simplicial_arrangement/simplicial_arrangement.h>

#include <spdlog/spdlog.h>

#include <catch2/catch.hpp>
#include <random>
#include <vector>

TEST_CASE("arrangement_scaling", "[ar][scaling][.benchmark]")
{
    using namespace simplicial_arrangement;
    std::mt19937 gen(7);
    std::uniform_int_distribution<> distrib(-(1 << 12), 1 << 12);

    for (size_t n = 10; n <= 100; n += 10) {
        BENCHMARK_ADVANCED(fmt::format("3D ar {} planes", n))
        (Catch::Benchmark::Chronometer meter)
        {
            disable_lookup_table();
            std::vector<Plane<double, 3>> planes;
            planes.reserve(n);
            for (size_t i = 0; i < n; i++) {
                planes.push_back({static_cast<double>(distrib(gen)),
                    static_cast<double>(distrib(gen)),
                    static_cast<double>(distrib(gen)),
                    static_cast<double>(distrib(gen))});
            }

            meter.measure([&]() { return compute_arrangement(planes); });
        };
    }
}

TEST_CASE("material_interface_scaling", "[mi][scaling][.benchmark]")
{
    using namespace simplicial_arrangement;
    std::mt19937 gen(7);
    std::uniform_int_distribution<> distrib(-(1 << 12), 1 << 12);

    for (size_t n = 10; n <= 100; n += 10) {
        BENCHMARK_ADVANCED(fmt::format("3D mi {} planes", n))
        (Catch::Benchmark::Chronometer meter)
        {
            disable_lookup_table();
            std::vector<Material<double, 3>> planes;
            planes.reserve(n);
            for (size_t i = 0; i < n; i++) {
                planes.push_back({static_cast<double>(distrib(gen)),
                    static_cast<double>(distrib(gen)),
                    static_cast<double>(distrib(gen)),
                    static_cast<double>(distrib(gen))});
            }

            meter.measure([&]() { return compute_material_interface(planes); });
        };
    }
}
