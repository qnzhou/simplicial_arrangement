#include <array>
#include <fstream>
#include <iostream>
#include <vector>

#include <spdlog/spdlog.h>
#include <nlohmann/json.hpp>

#include <simplicial_arrangement/material_interface.h>

#ifdef LOOKUP_TABLE_GEN
namespace simplicial_arrangement {

extern std::array<std::array<bool, 3>, 4> vertex_signs; // { m0 > m1, m1 > m2, m0 > m2 } x 4.
extern std::array<bool, 6> edge_signs; // intersection of m0 and m1 > m2 for each edge.

} // namespace simplicial_arrangement

auto generate_lookup_table()
{
    using namespace simplicial_arrangement;
    constexpr size_t INVALID = std::numeric_limits<size_t>::max();
    std::vector<Material<Int, 3>> dummy_materials(3);

    // Special ordering so MaterialInterfaceBuilder will insert planes in order.
    dummy_materials[0] = {{1, 0, 0, 0}};
    dummy_materials[1] = {{3, 0, 0, 0}};
    dummy_materials[2] = {{2, 0, 0, 0}};

    const size_t num_configs_1 = 1 << 12;

    std::vector<MaterialInterface<3>> data;
    data.reserve(num_configs_1 * 2);
    std::vector<size_t> level_1_start_index;
    level_1_start_index.reserve(num_configs_1 + 1);

    auto infer_edge_sign = [](size_t v0, size_t v1) {
        if (v0 > v1) std::swap(v0, v1);

        // TODO: only check the entries that should be true.
        if (v0 == 0 && v1 == 0) return false;
        if (v0 == 0 && v1 == 1) return false;
        if (v0 == 0 && v1 == 2) return false;
        if (v0 == 0 && v1 == 5) return false;
        if (v0 == 0 && v1 == 6) return false;

        if (v0 == 1 && v1 == 1) return false;
        if (v0 == 1 && v1 == 2) return false;
        if (v0 == 1 && v1 == 5) return false;
        if (v0 == 1 && v1 == 7) return false;

        if (v0 == 2 && v1 == 2) return false;
        if (v0 == 2 && v1 == 6) return true;
        if (v0 == 2 && v1 == 7) return true;

        if (v0 == 5 && v1 == 5) return false;
        if (v0 == 5 && v1 == 6) return true;
        if (v0 == 5 && v1 == 7) return true;

        if (v0 == 6 && v1 == 6) return true;
        if (v0 == 6 && v1 == 7) return true;

        if (v0 == 7 && v1 == 7) return true;
        throw std::runtime_error("Invalid edge state!");
    };

    std::vector<size_t> ambiguous_edges;
    ambiguous_edges.reserve(6);

    for (size_t i = 0; i < num_configs_1; i++) {
        level_1_start_index.push_back(data.size());
        vertex_signs[0][0] = i & 1;
        vertex_signs[0][1] = i & 2;
        vertex_signs[0][2] = i & 4;
        vertex_signs[1][0] = i & 8;
        vertex_signs[1][1] = i & 16;
        vertex_signs[1][2] = i & 32;
        vertex_signs[2][0] = i & 64;
        vertex_signs[2][1] = i & 128;
        vertex_signs[2][2] = i & 256;
        vertex_signs[3][0] = i & 512;
        vertex_signs[3][1] = i & 1024;
        vertex_signs[3][2] = i & 2048;

        size_t v0 = i & 7;
        size_t v1 = (i >> 3) & 7;
        size_t v2 = (i >> 6) & 7;
        size_t v3 = (i >> 9) & 7;

        // Skip impossible cyclical_ordering.
        if (v0 == 3 || v0 == 4) continue;
        if (v1 == 3 || v1 == 4) continue;
        if (v2 == 3 || v2 == 4) continue;
        if (v3 == 3 || v3 == 4) continue;

        ambiguous_edges.clear();
        if (v0 + v1 == 7)
            ambiguous_edges.push_back(0);
        else
            edge_signs[0] = infer_edge_sign(v0, v1);
        if (v0 + v2 == 7)
            ambiguous_edges.push_back(1);
        else
            edge_signs[1] = infer_edge_sign(v0, v2);
        if (v0 + v3 == 7)
            ambiguous_edges.push_back(2);
        else
            edge_signs[2] = infer_edge_sign(v0, v3);
        if (v1 + v2 == 7)
            ambiguous_edges.push_back(3);
        else
            edge_signs[3] = infer_edge_sign(v1, v2);
        if (v1 + v3 == 7)
            ambiguous_edges.push_back(4);
        else
            edge_signs[4] = infer_edge_sign(v1, v3);
        if (v2 + v3 == 7)
            ambiguous_edges.push_back(5);
        else
            edge_signs[5] = infer_edge_sign(v2, v3);

        if (ambiguous_edges.empty()) {
            data.push_back(compute_material_interface(dummy_materials));
        } else if (ambiguous_edges.size() == 1) {
            edge_signs[ambiguous_edges[0]] = false;
            data.push_back(compute_material_interface(dummy_materials));
            edge_signs[ambiguous_edges[0]] = true;
            data.push_back(compute_material_interface(dummy_materials));
        } else if (ambiguous_edges.size() == 2) {
            size_t num_cases = 1 << 2;
            for (size_t j = 0; j < num_cases; j++) {
                edge_signs[ambiguous_edges[0]] = j & 1;
                edge_signs[ambiguous_edges[1]] = j & 2;
                data.push_back(compute_material_interface(dummy_materials));
            }
        } else if (ambiguous_edges.size() == 3) {
            size_t num_cases = 1 << 3;
            for (size_t j = 0; j < num_cases; j++) {
                edge_signs[ambiguous_edges[0]] = j & 1;
                edge_signs[ambiguous_edges[1]] = j & 2;
                edge_signs[ambiguous_edges[2]] = j & 4;
                data.push_back(compute_material_interface(dummy_materials));
            }
        } else if (ambiguous_edges.size() == 4) {
            // It is impossible for ambiguous to form a loop. (Need proof.)
            if (ambiguous_edges == std::vector<size_t>({0, 1, 4, 5})) continue;
            if (ambiguous_edges == std::vector<size_t>({0, 2, 3, 5})) continue;
            if (ambiguous_edges == std::vector<size_t>({1, 2, 3, 4})) continue;
            if (ambiguous_edges == std::vector<size_t>({1, 2, 3, 4})) continue;
            size_t num_cases = 1 << 4;
            for (size_t j = 0; j < num_cases; j++) {
                edge_signs[ambiguous_edges[0]] = j & 1;
                edge_signs[ambiguous_edges[1]] = j & 2;
                edge_signs[ambiguous_edges[2]] = j & 4;
                edge_signs[ambiguous_edges[3]] = j & 8;
                data.push_back(compute_material_interface(dummy_materials));
            }
        } else {
            std::cout << v0 << ", " << v1 << ", " << v2 << ", " << v3 << " with ";
            std::cout << ambiguous_edges.size() << " ambiguous edges." << std::endl;
        }
    }
    level_1_start_index.push_back(data.size());
    return std::make_tuple(data, level_1_start_index);
}

void serialize(const std::vector<simplicial_arrangement::MaterialInterface<3>>& data,
    const std::vector<size_t>& level_1_start_index)
{
    nlohmann::json json;
    json = {{"data", nlohmann::json::array()}, {"start_index", level_1_start_index}};

    auto serialize_mi = [](const auto& mi) {
        nlohmann::json j;
        j.push_back(mi.vertices);
        j.push_back(nlohmann::json::array());
        j.push_back(nlohmann::json::array());

        for (const auto& f : mi.faces) {
            j[1].push_back({f.vertices, f.positive_material_label, f.negative_material_label});
        }

        for (const auto& c : mi.cells) {
            j[2].push_back({c.faces, c.material_label});
        }
        return j;
    };

    for (const auto& mi : data) {
        json["data"].push_back(serialize_mi(mi));
    }

    const std::string filename("material_interface_lookup_table.msgpack");
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
    std::array<size_t, 3> num_cells = {0, 0, 0};
    for (const auto& mi : data) {
        num_cells.at(mi.cells.size() - 1)++;
    }
    std::cout << "1 cell: " << num_cells[0] << std::endl;
    std::cout << "2 cell: " << num_cells[1] << std::endl;
    std::cout << "3 cell: " << num_cells[2] << std::endl;
    size_t num_level_1_entries = 0;
    size_t num_level_2_entries = 0;
    std::array<size_t, 5> ambiguous_cases = {0, 0, 0, 0};
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
            } else if (diff == 16) {
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
