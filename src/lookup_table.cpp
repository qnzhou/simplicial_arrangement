#include <simplicial_arrangement/lookup_table.h>
#include <simplicial_arrangement/material_interface.h>
#include <simplicial_arrangement/simplicial_arrangement.h>
#include "common.h"

#include <implicit_predicates/implicit_predicates.h>
#include <nlohmann/json.hpp>

#include <bitset>
#include <fstream>
#include <iostream>

namespace simplicial_arrangement {

// Global variables.
std::vector<Arrangement<3>> ar_data;
std::vector<size_t> ar_indices;
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
bool load_material_interface_lookup_table(const std::string& filename);
bool load_simplicial_arrangement_lookup_table(const std::string& filename);

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
        std::string simplicial_arrangement_lookup_table_file =
            lookup_table_path + "/simplicial_arrangement_lookup_table.msgpack";
        if (!load_simplicial_arrangement_lookup_table(simplicial_arrangement_lookup_table_file))
            return false;
    }

    if (table_type & MATERIAL_INTERFACE) {
        std::string material_interface_lookup_table_file =
            lookup_table_path + "/material_interface_lookup_table.msgpack";
        if (!load_material_interface_lookup_table(material_interface_lookup_table_file))
            return false;
    }

    use_lookup_table = true;
    return true;
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

bool load_simplicial_arrangement_lookup_table(const std::string& filename)
{
    if (!ar_data.empty() && !ar_indices.empty()) return true;
    auto t0 = std::chrono::high_resolution_clock::now();
    std::ifstream fin(filename.c_str(), std::ios::in | std::ios::binary);
    if (!fin) {
        std::cout << "Simplicial arrangement lookup table file not exist!" << std::endl;
        return false;
    }
    std::vector<char> msgpack(std::istreambuf_iterator<char>(fin), {});
    nlohmann::json json = nlohmann::json::from_msgpack(msgpack);
    fin.close();

    const auto deserialize_ar = [](const nlohmann::json& entry) {
        Arrangement<3> ar;
        ar.vertices = entry[0].get<std::vector<Point<3>>>();
        assert(ar.vertices.size() > 0);

        ar.faces.reserve(entry[1].size());
        for (const auto& face_entry : entry[1]) {
            Arrangement<3>::Face f;
            f.vertices = face_entry[0].get<std::vector<size_t>>();
            f.supporting_plane = face_entry[1].get<size_t>();
            f.positive_cell = face_entry[2].get<size_t>();
            f.negative_cell = face_entry[3].get<size_t>();
            ar.faces.push_back(std::move(f));
        }
        assert(ar.faces.size() > 0);

        ar.cells.reserve(entry[2].size());
        for (const auto& cell_entry : entry[2]) {
            Arrangement<3>::Cell c;
            c.faces = cell_entry.get<std::vector<size_t>>();
            ar.cells.push_back(std::move(c));
        }
        assert(ar.cells.size() > 0);

        return ar;
    };

    ar_indices = json["start_index"].get<std::vector<size_t>>();
    ar_data.reserve(json["data"].size());
    for (const auto& entry : json["data"]) {
        ar_data.push_back(deserialize_ar(entry));
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    logger().info("Simplicial arrangement lookup table load time: {}s",
        std::chrono::duration<double>(t1 - t0).count());
    return true;
}

namespace {

template <typename Scalar>
size_t ar_compute_outer_index_impl(const Plane<Scalar, 3>& p0)
{
    // Plane must not intersect tet at vertices.
    if (p0[0] == 0 || p0[1] == 0 || p0[2] == 0 || p0[3] == 0) return INVALID;

    // To reuse the 2 plane lookup table, we assume the second plane is negative
    // on all tet vertices.
    size_t index = 0;
    if (p0[0] > 0) index |= 1;
    if (p0[1] > 0) index |= 4;
    if (p0[2] > 0) index |= 16;
    if (p0[3] > 0) index |= 64;

    return index;
}

template <typename Scalar>
size_t ar_compute_outer_index_impl(const Plane<Scalar, 3>& p0, const Plane<Scalar, 3>& p1)
{
    // Plane must not intersect tet at vertices.
    if (p0[0] == 0 || p0[1] == 0 || p0[2] == 0 || p0[3] == 0) return INVALID;
    if (p1[0] == 0 || p1[1] == 0 || p1[2] == 0 || p1[3] == 0) return INVALID;

    size_t index = 0;
    if (p0[0] > 0) index |= 1;
    if (p1[0] > 0) index |= 2;
    if (p0[1] > 0) index |= 4;
    if (p1[1] > 0) index |= 8;
    if (p0[2] > 0) index |= 16;
    if (p1[2] > 0) index |= 32;
    if (p0[3] > 0) index |= 64;
    if (p1[3] > 0) index |= 128;

    return index;
}

template <typename Scalar>
size_t ar_compute_inner_index_impl(
    size_t outer_index, const Plane<Scalar, 3>& p0, const Plane<Scalar, 3>& p1)
{
    std::bitset<2> v0 = outer_index & 3;
    std::bitset<2> v1 = (outer_index >> 2) & 3;
    std::bitset<2> v2 = (outer_index >> 4) & 3;
    std::bitset<2> v3 = (outer_index >> 6) & 3;

    size_t index = 0;
    size_t edge_count = 0;
    std::array<Scalar, 2> pp0, pp1;

    auto add_edge = [&](size_t i, size_t j) -> bool {
        pp0 = {p0[i], p0[j]};
        pp1 = {p1[i], p1[j]};
        const auto s = implicit_predicates::orient1d(pp0.data(), pp1.data());
        if (s == implicit_predicates::ZERO || s == implicit_predicates::INVALID) return false;

        if (s > 0) index |= (1 << edge_count);
        edge_count++;
        return true;
    };

    if ((v0 ^ v1).all())
        if (!add_edge(0, 1)) return INVALID;
    if ((v0 ^ v2).all())
        if (!add_edge(0, 2)) return INVALID;
    if ((v0 ^ v3).all())
        if (!add_edge(0, 3)) return INVALID;
    if ((v1 ^ v2).all())
        if (!add_edge(1, 2)) return INVALID;
    if ((v1 ^ v3).all())
        if (!add_edge(1, 3)) return INVALID;
    if ((v2 ^ v3).all())
        if (!add_edge(2, 3)) return INVALID;

    if (edge_count == 4) {
        assert(index != 6 && index != 9); // Impossible cases.
        if (index < 6) {
        } else if (index < 9) {
            index -= 1; // Skipping invalid case with index 6.
        } else {
            index -= 2; // Skipping invalid case with index 6 and 9.
        }
    }

    return index;
}

template <typename Scalar>
size_t mi_compute_outer_index_impl(const Material<Scalar, 3>& m0, const Material<Scalar, 3>& m1)
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
    if (m0[0] == m1[0] || m1[0] == m2[0] || m0[0] == m2[0]) return INVALID;
    if (m0[1] == m1[1] || m1[1] == m2[1] || m0[1] == m2[1]) return INVALID;
    if (m0[2] == m1[2] || m1[2] == m2[2] || m0[2] == m2[2]) return INVALID;
    if (m0[3] == m1[3] || m1[3] == m2[3] || m0[3] == m2[3]) return INVALID;

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

size_t ar_compute_outer_index(const Plane<Int, 3>& p0)
{
    return ar_compute_outer_index_impl(p0);
}
size_t ar_compute_outer_index(const Plane<double, 3>& p0)
{
    return ar_compute_outer_index_impl(p0);
}
size_t ar_compute_outer_index(const Plane<Int, 3>& p0, const Plane<Int, 3>& p1)
{
    return ar_compute_outer_index_impl(p0, p1);
}
size_t ar_compute_outer_index(const Plane<double, 3>& p0, const Plane<double, 3>& p1)
{
    return ar_compute_outer_index_impl(p0, p1);
}
size_t ar_compute_inner_index(size_t outer_index, const Plane<Int, 3>& p0, const Plane<Int, 3>& p1)
{
    return ar_compute_inner_index_impl(outer_index, p0, p1);
}
size_t ar_compute_inner_index(
    size_t outer_index, const Plane<double, 3>& p0, const Plane<double, 3>& p1)
{
    return ar_compute_inner_index_impl(outer_index, p0, p1);
}

size_t mi_compute_outer_index(const Material<Int, 3>& m0, const Material<Int, 3>& m1)
{
    return mi_compute_outer_index_impl(m0, m1);
}
size_t mi_compute_outer_index(const Material<double, 3>& m0, const Material<double, 3>& m1)
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
