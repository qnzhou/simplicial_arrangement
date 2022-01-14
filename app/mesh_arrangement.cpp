//
// Created by Charles Du on 1/4/22.
//

#include <iostream>
#include <numeric> //std::partial_sum

#include "implicit_arrangement_util.h"
#include <igl/marching_tets.h>
#include "PyMesh/Arrangement.h"

#include "ScopedTimer.h"
#include <CLI/CLI.hpp>


struct IGL_Mesh {
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> vertices;
    Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor> faces;
};

void merge_meshes(const std::vector<IGL_Mesh>& meshes,
    // output
    IGL_Mesh &merged_mesh,
    Eigen::VectorXi &face_to_mesh)
{
    size_t num_meshes = meshes.size();
    std::vector<size_t> num_vertices(num_meshes + 1, 0);
    std::vector<size_t> num_faces(num_meshes + 1, 0);

    for (size_t i = 0; i < num_meshes; i++) {
        num_vertices[i + 1] = meshes[i].vertices.rows();
        num_faces[i + 1] = meshes[i].faces.rows();
    }

    std::partial_sum(num_vertices.begin(), num_vertices.end(),
        num_vertices.begin());
    std::partial_sum(num_faces.begin(), num_faces.end(), num_faces.begin());

    merged_mesh.vertices.resize(num_vertices.back(), 3);
    merged_mesh.faces.resize(num_faces.back(), 3);

    // face id to mesh id
    face_to_mesh = Eigen::VectorXi::Ones(num_faces.back());


    for (size_t i = 0; i < num_meshes; i++) {
        merged_mesh.vertices.block(num_vertices[i], 0,
            num_vertices[i + 1] - num_vertices[i], 3) =
            meshes[i].vertices;
        merged_mesh.faces.block(num_faces[i], 0,
            num_faces[i + 1] - num_faces[i], 3) =
            meshes[i].faces.array() + num_vertices[i];
        //
        face_to_mesh.segment(num_faces[i],num_faces[i+1]-num_faces[i]) *= i;
    }

}

IGL_Mesh generate_cube(const std::array<double,3> & bbox_min, const std::array<double,3> & bbox_max) {
    IGL_Mesh cube;
    cube.vertices.resize(8, 3);
    cube.vertices <<
        bbox_min[0], bbox_min[1], bbox_max[2],
        bbox_min[0], bbox_min[1], bbox_min[2],
        bbox_max[0], bbox_min[1], bbox_min[2],
        bbox_max[0], bbox_min[1], bbox_max[2],
        bbox_min[0], bbox_max[1], bbox_max[2],
        bbox_max[0], bbox_max[1], bbox_max[2],
        bbox_max[0], bbox_max[1], bbox_min[2],
        bbox_min[0], bbox_max[1], bbox_min[2];

    cube.faces.resize(12, 3);
    cube.faces << 0, 1, 2, 0, 2, 3, 4, 5, 6, 4, 6, 7, 0, 3, 5, 0, 5, 4, 1, 7, 6,
        1, 6, 2, 2, 6, 5, 2, 5, 3, 0, 4, 7, 0, 7, 1;

    return cube;
}



int main(int argc, const char* argv[]) {
    struct {
        std::string config_file;
        bool timing_only = false;
    } args;
    CLI::App app{"Mesh Arrangement Command Line"};
    app.add_option("config_file", args.config_file, "Configuration file")
        ->required();
    app.add_option("-T,--timing-only", args.timing_only, "Record timing without output result" );
    CLI11_PARSE(app, argc, argv);

    // record timings
    std::vector<std::string> timing_labels;
    std::vector<double> timings;

    // parse configure file
    std::string tet_mesh_file;
    std::string sphere_file;
    std::string output_dir;
    std::array<double,3> bbox_min, bbox_max;
    bool use_bbox = true;
    bool b_place_holder;
    parse_config_file(args.config_file, tet_mesh_file, sphere_file, output_dir,
        b_place_holder, b_place_holder,
        use_bbox, bbox_min, bbox_max);

    // load tet mesh
    std::vector<std::array<double, 3>> pts;
    std::vector<std::array<size_t, 4>> tets;
    load_tet_mesh(tet_mesh_file, pts, tets);
    size_t n_pts = pts.size();
    size_t n_tets = tets.size();
    std::cout << "tet mesh: " << pts.size() << " verts, " << tets.size() << " tets." << std::endl;


    // load implicit function values, or evaluate
    std::vector<Sphere> spheres;
    load_spheres(sphere_file, spheres);
    size_t n_func = spheres.size();

    std::vector<Eigen::VectorXd> funcVals;
    {
        timing_labels.emplace_back("func values");
        ScopedTimer<> timer("func values");
        funcVals.resize(n_func);
        for (size_t i = 0; i < n_func; i++) {
            funcVals[i].resize(pts.size());
            for (size_t j = 0; j < pts.size(); j++) {
                funcVals[i][j] = sphere_function(spheres[i].first, spheres[i].second, pts[j]);
            }
        }
        timings.push_back(timer.toc());
    }


    // convert tet mesh and function values to Eigen matrix
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> TV;
    TV.resize(pts.size(), Eigen::NoChange);
    for (size_t i = 0; i < pts.size(); ++i) {
        TV(i,0) = pts[i][0];
        TV(i,1) = pts[i][1];
        TV(i,2) = pts[i][2];
    }
    Eigen::Matrix<int, Eigen::Dynamic, 4, Eigen::RowMajor> TT;
    TT.resize(tets.size(), Eigen::NoChange);
    for (size_t i = 0; i < tets.size(); ++i) {
        TT(i, 0) = tets[i][0];
        TT(i, 1) = tets[i][1];
        TT(i, 2) = tets[i][2];
        TT(i, 3) = tets[i][3];
    }

    // marching tet
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>> SV;
    std::vector<Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>> SF;
    {
        timing_labels.emplace_back("marching-tets");
        ScopedTimer<> timer("marching-tets");
        SV.resize(n_func);
        SF.resize(n_func);
        for (size_t i = 0; i < n_func; ++i) {
            igl::marching_tets(TV, TT, funcVals[i], 0, SV[i], SF[i]);
        }
        timings.push_back(timer.toc());
    }

    // test: export iso-meshes
    if (!args.timing_only)
    {
        // convert marching tet results back to c++ vector
        std::vector<std::vector<std::array<double, 3>>> iso_pts_list(n_func);
        for (size_t i = 0; i < n_func; ++i) {
            const auto& sv = SV[i];
            auto& iso_verts = iso_pts_list[i];
            iso_verts.resize(sv.rows());
            for (size_t j = 0; j < iso_verts.size(); ++j) {
                iso_verts[j][0] = sv(j, 0);
                iso_verts[j][1] = sv(j, 1);
                iso_verts[j][2] = sv(j, 2);
            }
        }
        std::vector<std::vector<std::array<size_t, 3>>> iso_faces_list(n_func);
        for (size_t i = 0; i < n_func; ++i) {
            const auto& sf = SF[i];
            auto& iso_faces = iso_faces_list[i];
            iso_faces.resize(sf.rows());
            for (size_t j = 0; j < iso_faces.size(); ++j) {
                iso_faces[j][0] = sf(j, 0);
                iso_faces[j][1] = sf(j, 1);
                iso_faces[j][2] = sf(j, 2);
            }
        }
        // export iso-mesh
        save_tri_mesh_list(output_dir + "/marching_tet_iso_mesh.json",
            iso_pts_list,
            iso_faces_list);
    }

    // merge meshes
    std::vector<IGL_Mesh> iso_meshes;
    IGL_Mesh merged_mesh;
    Eigen::VectorXi face_to_mesh;
    {
        timing_labels.emplace_back("merge meshes");
        ScopedTimer<> timer("merge meshes");
        iso_meshes.resize(n_func);
        for (size_t i = 0; i < n_func; ++i) {
            iso_meshes[i].vertices = SV[i];
            iso_meshes[i].faces = SF[i];
        }
        if (use_bbox) {
            iso_meshes.push_back(generate_cube(bbox_min, bbox_max));
        }
        merge_meshes(iso_meshes, merged_mesh, face_to_mesh);
        timings.push_back(timer.toc());
    }

    // compute arrangement
    std::shared_ptr<PyMesh::Arrangement> engine;
    {
        ScopedTimer<> timer("mesh arrangement");
        engine = PyMesh::Arrangement::create_fast_arrangement(
            merged_mesh.vertices, merged_mesh.faces, face_to_mesh);
//        engine = PyMesh::Arrangement::create_mesh_arrangement(
//            merged_mesh.vertices, merged_mesh.faces, face_to_mesh);
//        engine->run();
        engine->run_with_timer(timing_labels, timings);
        timings.push_back(timer.toc());
    }
    timing_labels.emplace_back("arrangement(other)");
    double resolve_time = timings[timings.size()-3];
    double extract_time = timings[timings.size()-2];
    timings.back() = timings.back() - extract_time - resolve_time;
    std::cout << "Arrangement: resolving self-intersection: " << resolve_time << " s" <<  std::endl;
    std::cout << "Arrangement: extracting arrangement: " << extract_time << " s" << std::endl;


    // test: export arrangement cells
    size_t num_cells = engine->get_num_cells();
    std::cout << "num_cells = " << num_cells << std::endl;
    if (!args.timing_only)
    {
        // convert mesh arrangement cells to c++ vector
        std::vector<std::vector<std::array<double, 3>>> iso_pts_list(num_cells);
        for (size_t i = 0; i < num_cells; ++i) {
            const auto& sv = engine->get_vertices();
            auto& iso_verts = iso_pts_list[i];
            iso_verts.resize(sv.rows());
            for (size_t j = 0; j < iso_verts.size(); ++j) {
                iso_verts[j][0] = sv(j, 0);
                iso_verts[j][1] = sv(j, 1);
                iso_verts[j][2] = sv(j, 2);
            }
        }
        std::vector<std::vector<std::array<size_t, 3>>> iso_faces_list(num_cells);
        for (size_t i = 0; i < num_cells; ++i) {
            const auto& sf = engine->get_cell_faces(i);
            auto& iso_faces = iso_faces_list[i];
            iso_faces.resize(sf.rows());
            for (size_t j = 0; j < iso_faces.size(); ++j) {
                iso_faces[j][0] = sf(j, 0);
                iso_faces[j][1] = sf(j, 1);
                iso_faces[j][2] = sf(j, 2);
            }
        }
        // export iso-mesh
        save_tri_mesh_list(
                output_dir + "/mesh_arr_cells.json", iso_pts_list, iso_faces_list);
    }

    // export timings
    save_timings(output_dir + "/mesh_timings.json",
        timing_labels, timings);


    return 0;
}