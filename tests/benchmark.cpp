#include <simplicial_arrangement/SimplicialArrangement.h>

#include <catch2/catch.hpp>
#include <random>
#include <vector>

TEST_CASE("benchmark", "[arrangement][.benchmark]")
{
    using namespace simplicial_arrangement;

    constexpr size_t N = 5e2;
    std::vector<Plane<Int, 2>> int_data;
    std::vector<Plane<double, 2>> double_data;
    int_data.reserve(N);
    double_data.reserve(N);

    std::mt19937 gen(7);
    std::uniform_int_distribution<> distrib(-(1 << 12), 1 << 12);

    for (size_t i = 0; i < N; i++) {
        Int v0, v1, v2;
        v0 = distrib(gen);
        v1 = distrib(gen);
        v2 = distrib(gen);
        int_data.push_back({v0, v1, v2});
        double_data.push_back(
            {static_cast<double>(v0), static_cast<double>(v1), static_cast<double>(v2)});
    }

    std::function<size_t(const BSPNode<2>&)> count_num_cells;
    count_num_cells = [&](const BSPNode<2>& root) -> size_t {
        if (root.negative != nullptr && root.positive != nullptr) {
            return count_num_cells(*root.negative) + count_num_cells(*root.positive);
        } else {
            return 1;
        }
    };

    BENCHMARK_ADVANCED("2D arrangement with int")(Catch::Benchmark::Chronometer meter)
    {
        SimplicialArrangement<Int, 2> int_arrangement;
        int_arrangement.set_planes(int_data);
        meter.measure([&](){
                int_arrangement.initialize();
                });
        REQUIRE(count_num_cells(int_arrangement.get_root()) == 31310);
    };

    BENCHMARK_ADVANCED("2D arrangement with double")(Catch::Benchmark::Chronometer meter)
    {
        SimplicialArrangement<double, 2> double_arrangement;
        double_arrangement.set_planes(double_data);
        meter.measure([&](){
                double_arrangement.initialize();
                });
        REQUIRE(count_num_cells(double_arrangement.get_root()) == 31310);
    };
}
