#include "look_up_table.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>

namespace simplicial_arrangement {

std::unique_ptr<std::vector<Arrangement<3>>> one_func_lookup_table;
std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> to_check_edge_table;
std::unique_ptr<std::vector<std::vector<Arrangement<3>>>> two_func_lookup_table;
bool use_lookup_table;

bool load_lookup_table() {
    use_lookup_table = false;
#ifndef LOOKUP_TABLE_PATH
    return false;
#endif // !LOOKUP_TABLE_PATH
    std::string lookup_table_path = LOOKUP_TABLE_PATH;

    std::string one_func_table_file = lookup_table_path + "/1_func_lookup_table.json";
    if (!load_one_func_lookup_table(one_func_table_file)) return false;

    std::string to_check_edge_table_file = lookup_table_path + "/to_check_edge_table.json";
    if (!load_to_check_edge_table(to_check_edge_table_file)) return false;

    std::string two_func_table_file = lookup_table_path + "/2_func_lookup_table.json";
    if(!load_two_func_lookup_table(two_func_table_file)) return false;

    use_lookup_table = true;
    return true;
}


bool load_one_func_lookup_table(const std::string& filename)
{
    using json = nlohmann::json;
    std::ifstream fin(filename.c_str());
    if (!fin) {
        std::cout << "One function lookup table file not exist!" << std::endl;
        return false;
    }
    json input;
    fin >> input;
    fin.close();

    // populate the lookup table
    size_t num_entry = 16; // 2^4
    one_func_lookup_table = std::make_unique<std::vector<Arrangement<3>>>(num_entry);
    for (size_t i = 0; i < num_entry; i++) {
        auto &arrangement = (*one_func_lookup_table)[i];
        load_arrangement(arrangement, input[i]);
    }
    return true;
}

bool load_to_check_edge_table(const std::string& filename)
{
    using json = nlohmann::json;
    std::ifstream fin(filename.c_str());
    if (!fin) {
        std::cout << "Two function to-check edge lookup table file not exist!" << std::endl;
        return false;
    }
    json input;
    fin >> input;
    fin.close();

    // populate the lookup table
    size_t num_entry = 256; // 2^(4+4)
    to_check_edge_table =
        std::make_unique<std::vector<std::vector<std::pair<int, int>>>>(num_entry);
    for (size_t i = 0; i < num_entry; i++) {
        size_t n_edge = input[i].size();
        if (n_edge > 0) {
            auto& edges = (*to_check_edge_table)[i];
            edges.resize(n_edge);
            for (size_t j = 0; j < n_edge; j++) {
                edges[j].first = input[i][j][0];
                edges[j].second = input[i][j][1];
            }
        }
    }
    return true;
}

bool load_two_func_lookup_table(const std::string& filename)
{
    using json = nlohmann::json;
    std::ifstream fin(filename.c_str());
    if (!fin) {
        std::cout << "Two function lookup table file not exist!" << std::endl;
        return false;
    }
    json input;
    fin >> input;
    fin.close();

    // populate the lookup table
    size_t num_entry = 256; // 2^(4+4)
    two_func_lookup_table =
        std::make_unique<std::vector<std::vector<Arrangement<3>>>>(num_entry);
    for (size_t i = 0; i < num_entry; i++) {
        size_t n_sub_entry = input[i].size();
        auto& arrangements = (*two_func_lookup_table)[i];
        arrangements.resize(n_sub_entry);
        for (size_t j = 0; j < n_sub_entry; j++) {
            load_arrangement(arrangements[j], input[i][j]);
        } 
    }
    return true;
}

void load_arrangement(Arrangement<3>& arrangement, const nlohmann::json& data) {
    //
    auto& vertices = arrangement.vertices;
    vertices.resize(data[0].size());
    for (size_t j = 0; j < vertices.size(); j++) {
        vertices[j] = data[0][j];
    }
    //
    auto& faces = arrangement.faces;
    faces.resize(data[1].size());
    for (size_t j = 0; j < faces.size(); j++) {
        auto& face = faces[j];
        load_vector(face.vertices, data[1][j]);
        face.supporting_plane = data[2][j];
        face.positive_cell = data[3][j] < 0 ? Arrangement<3>::None : data[3][j];
        face.negative_cell = data[4][j] < 0 ? Arrangement<3>::None : data[4][j];
    }
    //
    auto& cells = arrangement.cells;
    cells.resize(data[5].size());
    for (size_t j = 0; j < cells.size(); j++) {
        auto& cell = cells[j];
        load_vector(cell.faces, data[5][j]);
        load_vector(cell.face_orientations, data[6][j]);
        load_vector(cell.plane_orientations, data[7][j]);
    }
    //
    load_vector(arrangement.unique_plane_indices, data[8]);
    auto& unique_planes = arrangement.unique_planes;
    unique_planes.resize(data[9].size());
    for (size_t j = 0; j < unique_planes.size(); j++) {
        load_vector(unique_planes[j], data[9][j]);
    }
    load_vector(arrangement.unique_plane_orientations, data[10]);
    //done
}

}