#include <simplicial_arrangement/lookup_table.h>
#include <simplicial_arrangement/material_interface.h>
#include <simplicial_arrangement/simplicial_arrangement.h>
#include "common.h"

#include <implicit_predicates/implicit_predicates.h>

#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

namespace simplicial_arrangement {

namespace {

template <typename T>
void load_vector(std::vector<T>& vec, const nlohmann::json& data)
{
    vec.resize(data.size());
    for (size_t i = 0; i < vec.size(); i++) {
        vec[i] = data[i].get<T>();
    }
}

} // namespace

std::unique_ptr<std::vector<Arrangement<3>>> one_func_lookup_table;
std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> to_check_edge_table;
std::unique_ptr<std::vector<std::vector<Arrangement<3>>>> two_func_lookup_table;
std::vector<MaterialInterface<3>> mi_data;
std::vector<size_t> mi_indices;
bool use_lookup_table = false;

#ifdef LOOKUP_TABLE_GEN
// Lookup-table-generation-related variables.
// These global varaibles are for predefined `mi_cut_0_face` function.
std::array<std::array<bool, 3>, 4> vertex_signs;
std::array<bool, 6> edge_signs;
#endif

// Forward declarations.
bool load_one_func_lookup_table(const std::string& filename);
bool load_to_check_edge_table(const std::string& filename);
bool load_two_func_lookup_table(const std::string& filename);
void load_arrangement(Arrangement<3>& arrangement, const nlohmann::json& data);
bool load_material_interface_lookup_table(const std::string& filename);

void enable_lookup_table()
{
    use_lookup_table = true;
}
void disable_lookup_table()
{
    use_lookup_table = false;
}

bool load_lookup_table(LookupTableType table_type)
{
    use_lookup_table = false;
#ifndef LOOKUP_TABLE_PATH
    return false;
#endif // !LOOKUP_TABLE_PATH
    std::string lookup_table_path = LOOKUP_TABLE_PATH;

    if (table_type & ARRANGEMENT) {
        std::string one_func_table_file = lookup_table_path + "/1_func_lookup_table.json";
        if (!load_one_func_lookup_table(one_func_table_file)) return false;

        std::string to_check_edge_table_file = lookup_table_path + "/to_check_edge_table.json";
        if (!load_to_check_edge_table(to_check_edge_table_file)) return false;

        std::string two_func_table_file = lookup_table_path + "/2_func_lookup_table.json";
        if (!load_two_func_lookup_table(two_func_table_file)) return false;
    }

    if (table_type & MATERIAL_INTERFACE) {
        std::string material_interface_lookup_table_file =
            lookup_table_path + "/material_interface_lookup_table.msgpack";
        if (!load_material_interface_lookup_table(material_interface_lookup_table_file)) return false;
    }

    use_lookup_table = true;
    return true;
}


bool load_one_func_lookup_table(const std::string& filename)
{
    if (one_func_lookup_table != nullptr) return true;
    auto t0 = std::chrono::high_resolution_clock::now();
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
        auto& arrangement = (*one_func_lookup_table)[i];
        load_arrangement(arrangement, input[i]);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    logger().info("Arrangement 1 function lookup table load time: {}s",
            std::chrono::duration<double>(t1 - t0).count());
    return true;
}

bool load_to_check_edge_table(const std::string& filename)
{
    if (to_check_edge_table != nullptr) return true;
    auto t0 = std::chrono::high_resolution_clock::now();
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
    auto t1 = std::chrono::high_resolution_clock::now();
    logger().info("Arrangement check edge lookup table load time: {}s",
            std::chrono::duration<double>(t1 - t0).count());
    return true;
}

bool load_two_func_lookup_table(const std::string& filename)
{
    if (two_func_lookup_table != nullptr) return true;
    auto t0 = std::chrono::high_resolution_clock::now();
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
    two_func_lookup_table = std::make_unique<std::vector<std::vector<Arrangement<3>>>>(num_entry);
    for (size_t i = 0; i < num_entry; i++) {
        size_t n_sub_entry = input[i].size();
        auto& arrangements = (*two_func_lookup_table)[i];
        arrangements.resize(n_sub_entry);
        for (size_t j = 0; j < n_sub_entry; j++) {
            if (input[i][j].size() > 0) {
                load_arrangement(arrangements[j], input[i][j]);
            }
        }
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    logger().info("Arrangement 2 function lookup table load time: {}s",
            std::chrono::duration<double>(t1 - t0).count());
    return true;
}

void load_arrangement(Arrangement<3>& arrangement, const nlohmann::json& data)
{
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
        face.positive_cell = data[3][j] < 0 ? Arrangement<3>::None : data[3][j].get<size_t>();
        face.negative_cell = data[4][j] < 0 ? Arrangement<3>::None : data[4][j].get<size_t>();
    }
    //
    auto& cells = arrangement.cells;
    cells.resize(data[5].size());
    for (size_t j = 0; j < cells.size(); j++) {
        auto& cell = cells[j];
        load_vector(cell.faces, data[5][j]);
    }
    //
    load_vector(arrangement.unique_plane_indices, data[8]);
    auto& unique_planes = arrangement.unique_planes;
    unique_planes.resize(data[9].size());
    for (size_t j = 0; j < unique_planes.size(); j++) {
        load_vector(unique_planes[j], data[9][j]);
    }
    load_vector(arrangement.unique_plane_orientations, data[10]);
    // done
}

bool load_material_interface_lookup_table(const std::string& filename)
{
    if (!mi_data.empty() && !mi_indices.empty()) return true;
    auto t0 = std::chrono::high_resolution_clock::now();
    std::ifstream fin(filename.c_str(), std::ios::in | std::ios::binary);
    if (!fin) {
        std::cout << "Material interface lookup table file not exist!" << std::endl;
        return false;
    }
    std::vector<char> msgpack(std::istreambuf_iterator<char>(fin), {});
    nlohmann::json json = nlohmann::json::from_msgpack(msgpack);
    fin.close();

    const auto deserialize_mi = [](const nlohmann::json& entry) {
        MaterialInterface<3> mi;
        mi.vertices = entry[0].get<std::vector<Joint<3>>>();

        mi.faces.reserve(entry[1].size());
        for (const auto& face_entry : entry[1]) {
            MaterialInterface<3>::Face f;
            f.vertices = face_entry[0].get<std::vector<size_t>>();
            f.positive_material_label = face_entry[1].get<size_t>();
            f.negative_material_label = face_entry[2].get<size_t>();
            mi.faces.push_back(std::move(f));
        }

        mi.cells.reserve(entry[2].size());
        for (const auto& cell_entry : entry[2]) {
            MaterialInterface<3>::Cell c;
            c.faces = cell_entry[0].get<std::vector<size_t>>();
            c.material_label = cell_entry[1].get<size_t>();
            mi.cells.push_back(std::move(c));
        }

        return mi;
    };

    mi_indices = json["start_index"].get<std::vector<size_t>>();
    mi_data.reserve(json["data"].size());
    for (const auto& entry : json["data"]) {
        mi_data.push_back(deserialize_mi(entry));
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    logger().info("Material interface lookup table load time: {}s",
            std::chrono::duration<double>(t1 - t0).count());
    return true;
}

namespace {

template <typename Scalar>
size_t mi_compute_outer_index_impl(
    const Material<Scalar, 3>& m0, const Material<Scalar, 3>& m1)
{
    // Material ordering at tet vertices must be unique.
    if (m0[0] == m1[0]) return INVALID;
    if (m0[1] == m1[1]) return INVALID;
    if (m0[2] == m1[2]) return INVALID;
    if (m0[3] == m1[3]) return INVALID;

    size_t index = 3510; // = 2 | 4 | 16 | 32 | 128 | 256 | 1024 | 2048

    // To reuse the 3 material lookup table, we assume we have an imaginary
    // material, m2, which is smaller than both m0 and m1.

    if (m0[0] > m1[0]) index |= 1;

    if (m0[1] > m1[1]) index |= 8;

    if (m0[2] > m1[2]) index |= 64;

    if (m0[3] > m1[3]) index |= 512;

    return index;
}

template <typename Scalar>
size_t mi_compute_outer_index_impl(
    const Material<Scalar, 3>& m0, const Material<Scalar, 3>& m1, const Material<Scalar, 3>& m2)
{
    // Material ordering at tet vertices must be unique.
    if (m0[0] == m1[0]  || m1[0] == m2[0] || m0[0] == m2[0]) return INVALID;
    if (m0[1] == m1[1]  || m1[1] == m2[1] || m0[1] == m2[1]) return INVALID;
    if (m0[2] == m1[2]  || m1[2] == m2[2] || m0[2] == m2[2]) return INVALID;
    if (m0[3] == m1[3]  || m1[3] == m2[3] || m0[3] == m2[3]) return INVALID;

    size_t index = 0;

    if (m0[0] > m1[0]) index |= 1;
    if (m1[0] > m2[0]) index |= 2;
    if (m0[0] > m2[0]) index |= 4;

    if (m0[1] > m1[1]) index |= 8;
    if (m1[1] > m2[1]) index |= 16;
    if (m0[1] > m2[1]) index |= 32;

    if (m0[2] > m1[2]) index |= 64;
    if (m1[2] > m2[2]) index |= 128;
    if (m0[2] > m2[2]) index |= 256;

    if (m0[3] > m1[3]) index |= 512;
    if (m1[3] > m2[3]) index |= 1024;
    if (m0[3] > m2[3]) index |= 2048;

    return index;
}

template <typename Scalar>
size_t mi_compute_inner_index_impl(size_t outer_index,
    const Material<Scalar, 3>& m0,
    const Material<Scalar, 3>& m1,
    const Material<Scalar, 3>& m2)
{
    const size_t v0 = outer_index & 7;
    const size_t v1 = (outer_index >> 3) & 7;
    const size_t v2 = (outer_index >> 6) & 7;
    const size_t v3 = (outer_index >> 9) & 7;

    size_t index = 0;
    size_t edge_count = 0;
    std::array<Scalar, 2> mm0, mm1, mm2;

    auto add_edge = [&]() {
        const auto s = implicit_predicates::mi_orient1d(mm0.data(), mm1.data(), mm2.data());
        if (s == implicit_predicates::ZERO || s == implicit_predicates::INVALID) return false;
        if (s > 0) {
            index |= (1 << edge_count);
        }
        edge_count++;
        return true;
    };

    if (v0 + v1 == 7) {
        mm0 = {m0[0], m0[1]};
        mm1 = {m1[0], m1[1]};
        mm2 = {m2[0], m2[1]};
        if (!add_edge()) return INVALID;
    }
    if (v0 + v2 == 7) {
        mm0 = {m0[0], m0[2]};
        mm1 = {m1[0], m1[2]};
        mm2 = {m2[0], m2[2]};
        if (!add_edge()) return INVALID;
    }
    if (v0 + v3 == 7) {
        mm0 = {m0[0], m0[3]};
        mm1 = {m1[0], m1[3]};
        mm2 = {m2[0], m2[3]};
        if (!add_edge()) return INVALID;
    }
    if (v1 + v2 == 7) {
        mm0 = {m0[1], m0[2]};
        mm1 = {m1[1], m1[2]};
        mm2 = {m2[1], m2[2]};
        if (!add_edge()) return INVALID;
    }
    if (v1 + v3 == 7) {
        mm0 = {m0[1], m0[3]};
        mm1 = {m1[1], m1[3]};
        mm2 = {m2[1], m2[3]};
        if (!add_edge()) return INVALID;
    }
    if (v2 + v3 == 7) {
        mm0 = {m0[2], m0[3]};
        mm1 = {m1[2], m1[3]};
        mm2 = {m2[2], m2[3]};
        if (!add_edge()) return INVALID;
    }
    return index;
}

} // namespace

size_t mi_compute_outer_index(
    const Material<Int, 3>& m0, const Material<Int, 3>& m1)
{
    return mi_compute_outer_index_impl(m0, m1);
}
size_t mi_compute_outer_index(
    const Material<double, 3>& m0, const Material<double, 3>& m1)
{
    return mi_compute_outer_index_impl(m0, m1);
}
size_t mi_compute_outer_index(
    const Material<Int, 3>& m0, const Material<Int, 3>& m1, const Material<Int, 3>& m2)
{
    return mi_compute_outer_index_impl(m0, m1, m2);
}
size_t mi_compute_outer_index(
    const Material<double, 3>& m0, const Material<double, 3>& m1, const Material<double, 3>& m2)
{
    return mi_compute_outer_index_impl(m0, m1, m2);
}
size_t mi_compute_inner_index(size_t outer_index,
    const Material<Int, 3>& m0,
    const Material<Int, 3>& m1,
    const Material<Int, 3>& m2)
{
    return mi_compute_inner_index_impl(outer_index, m0, m1, m2);
}
size_t mi_compute_inner_index(size_t outer_index,
    const Material<double, 3>& m0,
    const Material<double, 3>& m1,
    const Material<double, 3>& m2)
{
    return mi_compute_inner_index_impl(outer_index, m0, m1, m2);
}

} // namespace simplicial_arrangement
