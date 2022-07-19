#include <simplicial_arrangement/lookup_table.h>
#include <simplicial_arrangement/material_interface.h>
#include <simplicial_arrangement/simplicial_arrangement.h>

#include "io.h"

#include <absl/container/flat_hash_map.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <CLI/CLI.hpp>
#include <nlohmann/json.hpp>

#include <array>
#include <fstream>
#include <string>

namespace {

template <typename T>
inline T det3(T a, T b, T c, T d, T e, T f, T g, T h, T i)
{
    return a * (e * i - f * h) - b * (d * i - g * f) + c * (d * h - e * g);
}

void save_edges(const std::string& filename,
    const std::vector<std::array<double, 4>>& points,
    const std::vector<std::array<size_t, 2>>& edges)
{
    wmtk::MshData data;
    data.add_edge_vertices(points.size(), [&](size_t i) {
        const auto& p = points[i];
        return std::array<double, 3>{p[1], p[2], p[3]};
    });
    data.add_edges(edges.size(), [&](size_t i) { return edges[i]; });
    data.save(filename);
}

template <typename F>
auto extract_triangles(const std::vector<F>& faces, const std::vector<bool>& selected)
{
    std::vector<std::array<size_t, 3>> triangles;
    std::vector<size_t> face_id;
    std::vector<size_t> supporting_planes;
    triangles.reserve(faces.size() * 2);
    face_id.reserve(faces.size() * 2);
    supporting_planes.reserve(faces.size() * 2);
    size_t fid = 0;
    for (const auto& f : faces) {
        if (selected[fid]) {
            const size_t s = f.vertices.size();
            for (size_t i = 1; i < s - 1; i++) {
                triangles.push_back({f.vertices[0], f.vertices[i], f.vertices[i + 1]});
                face_id.push_back(fid);
                supporting_planes.push_back(f.supporting_plane);
            }
        }
        fid++;
    }

    return std::make_tuple(triangles, face_id, supporting_planes);
}

template <typename F>
auto extract_triangles(const std::vector<F>& faces)
{
    std::vector<bool> selected(faces.size(), true);
    return extract_triangles(faces, selected);
}

template <typename F>
void save_faces(const std::string& filename,
    const std::vector<std::array<double, 4>>& points,
    const std::vector<F>& faces)
{
    wmtk::MshData data;
    data.add_face_vertices(points.size(), [&](size_t i) {
        const auto& p = points[i];
        return std::array<double, 3>{p[1], p[2], p[3]};
    });

    auto r = extract_triangles(faces);
    const auto& triangles = std::get<0>(r);
    const auto& face_id = std::get<1>(r);
    const auto& supporting_planes = std::get<2>(r);

    data.add_faces(triangles.size(), [&](size_t i) { return triangles[i]; });
    data.add_face_attribute<1>("id", [&](size_t i) { return face_id[i]; });
    data.add_face_attribute<1>("supporting_plane", [&](size_t i) { return supporting_planes[i]; });
    data.save(filename);
}

template <typename F, typename C>
auto extract_cell(
    const std::vector<std::array<double, 4>>& points, const std::vector<F>& faces, const C& cell)
{
    const size_t num_points = points.size();
    const size_t num_faces = faces.size();

    std::vector<bool> active_points(num_points, false);
    std::vector<bool> active_faces(num_faces, false);

    for (size_t fid : cell.faces) {
        active_faces[fid] = true;
        const auto& f = faces[fid];
        for (size_t vid : f.vertices) {
            active_points[vid] = true;
        }
    }

    auto r = extract_triangles(faces, active_faces);
    auto& triangles = std::get<0>(r);
    auto& face_id = std::get<1>(r);
    auto& supporting_planes = std::get<2>(r);

    std::vector<std::array<double, 3>> vertices;
    vertices.reserve(num_points);
    std::vector<size_t> vertex_map(num_points, num_points + 1);
    for (size_t i = 0; i < num_points; i++) {
        if (active_points[i]) {
            const auto& p = points[i];
            vertex_map[i] = vertices.size();
            vertices.push_back({p[1], p[2], p[3]});
        }
    }

    for (auto& t : triangles) {
        t[0] = vertex_map[t[0]];
        t[1] = vertex_map[t[1]];
        t[2] = vertex_map[t[2]];
    }
    return std::make_tuple(vertices, triangles, face_id, supporting_planes);
}

template <typename F, typename C>
void save_cells(const std::string& filename,
    const std::vector<std::array<double, 4>>& points,
    const std::vector<F>& faces,
    const std::vector<C>& cells)
{
    wmtk::MshData data;

    for (const auto& cell : cells) {
        auto r = extract_cell(points, faces, cell);
        const auto& vertices = std::get<0>(r);
        const auto& triangles = std::get<1>(r);
        data.add_face_vertices(vertices.size(), [&](size_t i) { return vertices[i]; });
        data.add_faces(triangles.size(), [&](size_t i) { return triangles[i]; });
    }
    data.save(filename);
}

} // namespace


spdlog::logger& logger()
{
    static auto default_logger = spdlog::stdout_color_mt("Simplex");
    return *default_logger;
}

auto load_config_file(const std::string& config_file)
    -> std::tuple<bool, std::vector<std::array<double, 4>>>
{
    logger().info("Loading {}", config_file);
    nlohmann::json data;
    std::ifstream fin(config_file.c_str());
    if (!fin.is_open()) {
        throw std::runtime_error("Config file not found!");
    }
    fin >> data;
    fin.close();

    bool is_mi = data["type"].get<std::string>() == "MI";
    size_t num_planes = data["planes"].size();
    std::vector<std::array<double, 4>> planes;
    planes.reserve(num_planes);

    for (const auto& plane_data : data["planes"]) {
        assert(plane_data.size() == 4);
        planes.push_back({plane_data[0].get<double>(),
            plane_data[1].get<double>(),
            plane_data[2].get<double>(),
            plane_data[3].get<double>()});
    }

    logger().info("type: {}", is_mi ? "MI" : "IA");
    for (const auto& p : planes) {
        logger().info("{}, {}, {}, {}", p[0], p[1], p[2], p[3]);
    }
    return {is_mi, planes};
}

auto compute_IA(const std::vector<std::array<double, 4>>& planes)
{
    auto r = simplicial_arrangement::compute_arrangement(planes);
    const size_t num_vertices = r.vertices.size();
    const size_t num_faces = r.faces.size();
    const size_t num_cells = r.cells.size();

    logger().info("[IA] #vertices: {}", num_vertices);
    logger().info("[IA]    #faces: {}", num_faces);
    logger().info("[IA]    #cells: {}", num_cells);

    const auto compute_barycentric_coords = [&](const simplicial_arrangement::Point<3>& p) {
        static std::array<std::array<double, 4>, 3> P;
        for (size_t i = 0; i < 3; i++) {
            if (p[i] < 4) {
                P[i] = {0, 0, 0, 0};
                P[i][p[i]] = 1;
            } else {
                assert(p[i] - 4 < planes.size());
                P[i] = planes[p[i] - 4];
            }
        }

        std::array<double, 4> b;
        // clang-format off
        b[0] = -::det3(
                P[0][1], P[0][2], P[0][3],
                P[1][1], P[1][2], P[1][3],
                P[2][1], P[2][2], P[2][3]);
        b[1] = ::det3(
                P[0][0], P[0][2], P[0][3],
                P[1][0], P[1][2], P[1][3],
                P[2][0], P[2][2], P[2][3]);
        b[2] = -::det3(
                P[0][0], P[0][1], P[0][3],
                P[1][0], P[1][1], P[1][3],
                P[2][0], P[2][1], P[2][3]);
        b[3] = ::det3(
                P[0][0], P[0][1], P[0][2],
                P[1][0], P[1][1], P[1][2],
                P[2][0], P[2][1], P[2][2]);
        double sum = b[0] + b[1] + b[2] + b[3];
        b[0] /= sum;
        b[1] /= sum;
        b[2] /= sum;
        b[3] /= sum;
        // clang-format on
        logger().info("barycentric coord: {}, {}, {}, {}", b[0], b[1], b[2], b[3]);
        return b;
    };

    const auto compute_nonmanifold_edges = [&]() {
        const auto& faces = r.faces;
        absl::flat_hash_map<std::array<size_t, 2>, size_t> edge_count;
        edge_count.reserve(faces.size() * 4);

        for (const auto& f : faces) {
            const size_t s = f.vertices.size();
            for (size_t i = 0; i < s; i++) {
                std::array<size_t, 2> e;
                e[0] = f.vertices[i];
                e[1] = f.vertices[(i + 1) % s];
                if (e[0] > e[1]) std::swap(e[0], e[1]);
                auto [itr, success] = edge_count.try_emplace(e, 1);
                if (!success) {
                    // edge exists.
                    itr->second++;
                }
            }
        }

        std::vector<std::array<size_t, 2>> nonmanifold_edges;
        nonmanifold_edges.reserve(edge_count.size());
        for (const auto& entry : edge_count) {
            // if (entry.second <= 2) continue;
            nonmanifold_edges.push_back(entry.first);
        }
        return nonmanifold_edges;
    };

    std::vector<std::array<double, 4>> vertices;
    vertices.reserve(num_vertices);
    for (const auto& p : r.vertices) {
        vertices.push_back(compute_barycentric_coords(p));
    }

    auto edges = compute_nonmanifold_edges();
    save_edges("ia_edges.msh", vertices, edges);
    save_faces("ia_faces.msh", vertices, r.faces);
    save_cells("ia_cells.msh", vertices, r.faces, r.cells);
}

int main(int argc, char** argv)
{
    struct
    {
        std::string config_file;
        std::string msh_file;
    } args;
    CLI::App app("IA/MI case study within a simplex.");
    app.add_option("config_file", args.config_file, "Configuration file")->required();
    CLI11_PARSE(app, argc, argv);

    using namespace simplicial_arrangement;
    spdlog::set_level(spdlog::level::info);

    auto [is_mi, planes] = load_config_file(args.config_file);

    if (is_mi) {
        // TODO
    } else {
        compute_IA(planes);
    }

    return 0;
}
