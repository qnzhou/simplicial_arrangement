#include <bitset>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include <spdlog/spdlog.h>
#include <nlohmann/json.hpp>

#include <simplicial_arrangement/simplicial_arrangement.h>

#ifdef LOOKUP_TABLE_GEN
namespace simplicial_arrangement {

/**
 * Signs of f0 and f1 on each tet vertex.
 *
 * @note Only the first 2 boolean entries are used at each vertex.
 */
extern std::array<std::array<bool, 3>, 4> vertex_signs;

/**
 * Zero crossing of f0 is on the positive side of f1 for each edge.
 */
extern std::array<bool, 6> edge_signs;

} // namespace simplicial_arrangement

auto generate_lookup_table()
{
    using namespace simplicial_arrangement;
    constexpr size_t INVALID = std::numeric_limits<size_t>::max();
    std::vector<Plane<Int, 3>> dummy_planes(2);

    dummy_planes[0] = {{1, 0, 0, 0}};
    dummy_planes[1] = {{2, 0, 0, 0}};

    const size_t num_configs_1 = 1 << 8;

    std::vector<Arrangement<3>> data;
    data.reserve(num_configs_1 * 2);
    std::vector<size_t> level_1_start_index;
    level_1_start_index.reserve(num_configs_1 + 1);

    auto infer_edge_sign = [](auto v0, auto v1) -> bool {
        if (!(v0[0] ^ v1[0])) return false; // f0 has no zero crossing. Don't care.
        if (!(v0[1] ^ v1[1])) return v0[1]; // f1 has no zero crossing. Use its sign.

        assert(false);
        return false; // All other cases are ambiguous.
    };

    std::vector<size_t> ambiguous_edges;
    ambiguous_edges.reserve(6);

    for (size_t i = 0; i < num_configs_1; i++) {
        level_1_start_index.push_back(data.size());
        vertex_signs[0][0] = i & 1;
        vertex_signs[0][1] = i & 2;
        vertex_signs[1][0] = i & 4;
        vertex_signs[1][1] = i & 8;
        vertex_signs[2][0] = i & 16;
        vertex_signs[2][1] = i & 32;
        vertex_signs[3][0] = i & 64;
        vertex_signs[3][1] = i & 128;

        std::bitset<2> v0 = i & 3;
        std::bitset<2> v1 = (i >> 2) & 3;
        std::bitset<2> v2 = (i >> 4) & 3;
        std::bitset<2> v3 = (i >> 6) & 3;

        ambiguous_edges.clear();
        if ((v0 ^ v1).all())
            ambiguous_edges.push_back(0);
        else
            edge_signs[0] = infer_edge_sign(v0, v1);
        if ((v0 ^ v2).all())
            ambiguous_edges.push_back(1);
        else
            edge_signs[1] = infer_edge_sign(v0, v2);
        if ((v0 ^ v3).all())
            ambiguous_edges.push_back(2);
        else
            edge_signs[2] = infer_edge_sign(v0, v3);
        if ((v1 ^ v2).all())
            ambiguous_edges.push_back(3);
        else
            edge_signs[3] = infer_edge_sign(v1, v2);
        if ((v1 ^ v3).all())
            ambiguous_edges.push_back(4);
        else
            edge_signs[4] = infer_edge_sign(v1, v3);
        if ((v2 ^ v3).all())
            ambiguous_edges.push_back(5);
        else
            edge_signs[5] = infer_edge_sign(v2, v3);

        if (ambiguous_edges.empty()) {
            data.push_back(compute_arrangement(dummy_planes));
        } else if (ambiguous_edges.size() == 1) {
            edge_signs[ambiguous_edges[0]] = false;
            data.push_back(compute_arrangement(dummy_planes));
            edge_signs[ambiguous_edges[0]] = true;
            data.push_back(compute_arrangement(dummy_planes));
        } else {
            const size_t num_ambiguous_edges = ambiguous_edges.size();
            const size_t num_cases = 1 << num_ambiguous_edges;
            assert(num_ambiguous_edges <= 4);

            for (size_t j = 0; j < num_cases; j++) {
                for (size_t k = 0; k < num_ambiguous_edges; k++) {
                    edge_signs[ambiguous_edges[k]] = j & (1 << k);
                }

                // 4 tet edges must form a loop that consists of 2 opposing edge
                // pairs.
                if (ambiguous_edges == std::vector<size_t>({0, 1, 4, 5})) {
                    if ((edge_signs[0] != edge_signs[1]) && (edge_signs[1] != edge_signs[5]) &&
                        (edge_signs[5] != edge_signs[4])) {
                        spdlog::info("1. i={} Skipping j={}", i, j);
                        continue;
                    }
                } else if (ambiguous_edges == std::vector<size_t>({0, 2, 3, 5})) {
                    if ((edge_signs[0] != edge_signs[2]) && (edge_signs[2] != edge_signs[5]) &&
                        (edge_signs[5] != edge_signs[3])) {
                        spdlog::info("2. i={} Skipping j={}", i, j);
                        continue;
                    }
                } else if (ambiguous_edges == std::vector<size_t>({1, 2, 3, 4})) {
                    if ((edge_signs[1] != edge_signs[2]) && (edge_signs[2] != edge_signs[4]) &&
                        (edge_signs[4] != edge_signs[3])) {
                        spdlog::info("3. i={} Skipping j={}", i, j);
                        continue;
                    }
                }

                data.push_back(compute_arrangement(dummy_planes));
            }
        }
    }
    level_1_start_index.push_back(data.size());
    return std::make_tuple(data, level_1_start_index);
}

void serialize(const std::vector<simplicial_arrangement::Arrangement<3>>& data,
    const std::vector<size_t>& level_1_start_index)
{
    nlohmann::json json;
    json = {{"data", nlohmann::json::array()}, {"start_index", level_1_start_index}};

    auto serialize_ia = [](const auto& ia) {
        nlohmann::json j;
        j.push_back(ia.vertices);
        j.push_back(nlohmann::json::array());
        j.push_back(nlohmann::json::array());

        for (const auto& f : ia.faces) {
            j[1].push_back({f.vertices, f.supporting_plane, f.positive_cell, f.negative_cell});
        }

        for (const auto& c : ia.cells) {
            j[2].push_back(c.faces);
        }
        return j;
    };

    for (const auto& ia : data) {
        json["data"].push_back(serialize_ia(ia));
    }

    const std::string filename("simplicial_arrangement_lookup_table.msgpack");
    std::vector<std::uint8_t> msgpack = nlohmann::json::to_msgpack(json);
    std::ofstream fl(filename.c_str(), std::ios::out | std::ios::binary);
    fl.write(reinterpret_cast<char*>(msgpack.data()), msgpack.size());
    std::cout << "Lookup table outputed to " << filename << std::endl;

    std::ofstream fout("tmp.json", std::ios::out);
    fout << std::setw(4) << json;
}

#endif

int main(int argc, char** argv)
{
#ifndef LOOKUP_TABLE_GEN
    std::cerr << "This command requires the -DLOOKUP_TABLE_GEN flag to be on!" << std::endl;
    return -1;
#else
    auto [data, start_index] = generate_lookup_table();
    std::cout << data.size() << std::endl;
    std::array<size_t, 4> num_cells = {0, 0, 0, 0};
    for (const auto& ia : data) {
        num_cells.at(ia.cells.size() - 1)++;
    }
    std::cout << "1 cell: " << num_cells[0] << std::endl;
    std::cout << "2 cell: " << num_cells[1] << std::endl;
    std::cout << "3 cell: " << num_cells[2] << std::endl;
    std::cout << "4 cell: " << num_cells[3] << std::endl;
    size_t num_level_1_entries = 0;
    size_t num_level_2_entries = 0;
    std::array<size_t, 5> ambiguous_cases = {0, 0, 0, 0, 0};
    for (size_t i = 1; i < start_index.size(); i++) {
        size_t diff = start_index[i] - start_index[i - 1];
        if (diff == 0) continue;
        if (diff == 1) {
            num_level_1_entries++;
            ambiguous_cases[0]++;
        } else {
            num_level_2_entries += diff;
            if (diff == 2) {
                ambiguous_cases[1]++;
            } else if (diff == 4) {
                ambiguous_cases[2]++;
            } else if (diff == 8) {
                ambiguous_cases[3]++;
            } else if (diff <= 16) {
                ambiguous_cases[4]++;
            } else {
                throw std::runtime_error("Bad case");
            }
        }
    }
    std::cout << "level 1 size: " << num_level_1_entries << std::endl;
    std::cout << "level 2 size: " << num_level_2_entries << std::endl;
    std::cout << "with 1 ambiguous edge: " << ambiguous_cases[1] << std::endl;
    std::cout << "with 2 ambiguous edge: " << ambiguous_cases[2] << std::endl;
    std::cout << "with 3 ambiguous edge: " << ambiguous_cases[3] << std::endl;
    std::cout << "with 4 ambiguous edge: " << ambiguous_cases[4] << std::endl;

    serialize(data, start_index);

    return 0;
#endif
}
