//
// Created by Charles Du on 1/4/22.
//

#include <iostream>
#include <numeric> //std::partial_sum

#include "implicit_arrangement_util.h"
#include <igl/marching_tets.h>
#include "PyMesh/Arrangement.h"

#include "ScopedTimer.h"

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



int main(int argc, const char* argv[]) {
//    std::string dataDir = "D:/research/simplicial_arrangement/data/";
    std::string dataDir = "/Users/charlesdu/Downloads/research/implicit_modeling/code/simplicial_arrangement/data/";
    std::string resolution = "100k";
    // load tet mesh
    std::string tet_mesh_file = dataDir + "tet_mesh_" + resolution + ".json";
    std::vector<std::array<double, 3>> pts;
    std::vector<std::array<size_t, 4>> tets;
    load_tet_mesh(tet_mesh_file, pts, tets);
    size_t n_pts = pts.size();
    size_t n_tets = tets.size();
    std::cout << "tet mesh: " << pts.size() << " verts, " << tets.size() << " tets." << std::endl;

    // load implicit function values, or evaluate
    size_t n_func = 4;
    std::vector<std::array<double, 3>> centers(n_func);
    centers[0] = {0, 0, 0};
    centers[1] = {0.5, 0, 0};
    centers[2] = {0, 0.5, 0};
    centers[3] = {0, 0, 0.5};
    double radius = 0.5;

    std::vector<Eigen::VectorXd> funcVals;
    {
        ScopedTimer<> timer("evaluate implicit functions at vertices");
        funcVals.resize(n_func);
        for (size_t i = 0; i < n_func; i++) {
            funcVals[i].resize(pts.size());
            for (size_t j = 0; j < pts.size(); j++) {
                funcVals[i][j] = sphere_function(centers[i], radius, pts[j]);
            }
        }
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
//    std::vector<Eigen::VectorXd> func_S(n_func);
//    for (size_t i = 0; i < n_func; ++i) {
//        auto& S = func_S[i];
//        S.resize(pts.size());
//        for (size_t j = 0; j < pts.size(); ++j) {
//            S(j) = funcVals[i][j];
//        }
//    }

    // marching tet
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>> SV;
    std::vector<Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>> SF;
    {
        ScopedTimer<> timer("marching tet");
        SV.resize(n_func);
        SF.resize(n_func);
        for (size_t i = 0; i < n_func; ++i) {
            igl::marching_tets(TV, TT, funcVals[i], 0, SV[i], SF[i]);
        }
    }

    // test: export iso-meshes
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
        save_tri_mesh_list(dataDir + "marching_tet_iso_mesh_" + resolution + ".json",
            iso_pts_list,
            iso_faces_list);
    }

    // merge meshes
    std::vector<IGL_Mesh> iso_meshes;
    IGL_Mesh merged_mesh;
    Eigen::VectorXi face_to_mesh;
    {
        ScopedTimer<> timer("merge meshes");
        iso_meshes.resize(n_func);
        for (size_t i = 0; i < n_func; ++i) {
            iso_meshes[i].vertices = SV[i];
            iso_meshes[i].faces = SF[i];
        }
        merge_meshes(iso_meshes, merged_mesh, face_to_mesh);
    }

    // compute arrangement
    std::shared_ptr<PyMesh::Arrangement> engine;
    {
        ScopedTimer<> timer("(fast) mesh arrangement");
        engine = PyMesh::Arrangement::create_fast_arrangement(
            merged_mesh.vertices, merged_mesh.faces, face_to_mesh);
//        engine = PyMesh::Arrangement::create_mesh_arrangement(
//            merged_mesh.vertices, merged_mesh.faces, face_to_mesh);
        engine->run();
    }

    // test: export arrangement cells
    {
        size_t num_cells = engine->get_num_cells();
        std::cout << "num_cells = " << num_cells << std::endl;
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
        save_tri_mesh_list(dataDir + "mesh_arr_cells_" + resolution + ".json",
            iso_pts_list,
            iso_faces_list);
    }


    return 0;
}