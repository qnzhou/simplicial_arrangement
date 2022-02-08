//
// Created by Charles Du on 1/20/22.
//
#include "ScopedTimer.h"
#include "implicit_arrangement_util.h"

#include <absl/container/flat_hash_set.h>
#include <CLI/CLI.hpp>
#include <Eigen/Core>

int main(int argc, const char* argv[])
{
    struct
    {
        std::string config_file;
        bool timing_only = false;
    } args;
    CLI::App app{"Dual Contouring (Material Interface) Command Line"};
    app.add_option("config_file", args.config_file, "Configuration file")->required();
    app.add_option("-T,--timing-only", args.timing_only, "Record timing without output result");
    CLI11_PARSE(app, argc, argv);

    // parse configure file
    std::string tet_mesh_file;
    std::string material_file;
    std::string output_dir;
    bool b_place_holder;
    parse_config_file_MI(args.config_file, tet_mesh_file, material_file, output_dir,
        b_place_holder, b_place_holder, b_place_holder);

    // record timings
    std::vector<std::string> timing_labels;
    std::vector<double> timings;

    // load tet mesh
    std::vector<std::array<double, 3>> pts;
    std::vector<std::array<size_t, 4>> tets;
    load_tet_mesh(tet_mesh_file, pts, tets);
    size_t n_tets = tets.size();
    size_t n_pts = pts.size();
    std::cout << "tet mesh: " << pts.size() << " verts, " << tets.size() << " tets." << std::endl;

    // load tet mesh and function values
//    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> funcVals;
//    load_tet_mesh_func(tet_mesh_file, pts, tets, funcVals);
//    size_t n_tets = tets.size();
//    size_t n_pts = pts.size();
//    std::cout << "tet mesh: " << pts.size() << " verts, " << tets.size() << " tets." << std::endl;



    // load implicit functions and compute function values at vertices
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> funcVals;
    if (load_functions(material_file, pts, funcVals)) {
        std::cout << "function loading finished." << std::endl;
    } else {
        std::cout << "function loading failed." << std::endl;
        return -2;
    }
    size_t n_func = funcVals.cols();

    // highest material at vertices
    std::vector<size_t> highest_material_at_vert;
    {
        timing_labels.emplace_back("highest func");
        ScopedTimer<> timer("highest func");
        highest_material_at_vert.reserve(n_pts);
        for (Eigen::Index i = 0; i < n_pts; i++) {
            double max = funcVals(i, 0);
            size_t max_id = 0;
            for (Eigen::Index j = 1; j < n_func; j++) {
                if (funcVals(i,j) > max) {
                    max = funcVals(i,j);
                    max_id = j;
                }
            }
            highest_material_at_vert.push_back(max_id);
        }
        timings.push_back(timer.toc());
    }

    // filter relevant materials in each tet
    // a tet is non-empty if there are some material interface in it
    std::vector<bool> has_intersection;
    size_t num_2func = 0;
    size_t num_3func = 0;
    size_t num_more_func = 0;
    {
        timing_labels.emplace_back("filter");
        ScopedTimer<> timer("filter");
        has_intersection.resize(n_tets, false);
        absl::flat_hash_set<size_t> materials;
        materials.reserve(3);
        for (size_t i = 0; i < n_tets; ++i) {
            const auto& tet = tets[i];
            // find high materials
            materials.clear();
            for (size_t j = 0; j < 4; ++j) {
                materials.insert(highest_material_at_vert[tet[j]]);
            }
            // if more than one high material, there is material interface
            if (materials.size() > 1) {  // has material interface
                has_intersection[i] = true;
                switch (materials.size()) {
                case 2: ++num_2func; break;
                case 3: ++num_3func; break;
                default: ++num_more_func; break;
                }
            }
        }
        timings.push_back(timer.toc());
    }

    // dual contouring
    std::vector<std::array<double, 3>> mesh_verts;
    std::vector<std::array<size_t, 3>> mesh_tris;
    {
        timing_labels.emplace_back("dual contouring");
        ScopedTimer<> timer("dual contouring");
        tet_dual_contouring_MI(num_2func,
            num_3func,
            num_more_func,
            pts,
            tets,
            funcVals,
            highest_material_at_vert,
            has_intersection,
            mesh_verts,
            mesh_tris);
        timings.push_back(timer.toc());
    }
    std::cout << "num mesh-vertices = " << mesh_verts.size() << std::endl;
    std::cout << "num mesh-triangles = " << mesh_tris.size() << std::endl;


    // convert mesh_tris to mesh_faces
    std::vector<PolygonFace> mesh_faces;
    {
        mesh_faces.reserve(mesh_tris.size());
        for (auto& mesh_tri : mesh_tris) {
            mesh_faces.emplace_back();
            auto& vert_indices = mesh_faces.back().vert_indices;
            vert_indices.insert(vert_indices.end(), mesh_tri.begin(), mesh_tri.end());
        }
    }

    //  compute mesh-edges and edge-face connectivity
    std::vector<Edge> mesh_edges;
    std::vector<std::vector<size_t>> edges_of_mesh_face;
    {
        timing_labels.emplace_back("edge-face connectivity");
        ScopedTimer<> timer("Edge-face connectivity");
        compute_mesh_edges(mesh_faces, edges_of_mesh_face, mesh_edges);
        timings.push_back(timer.toc());
    }
    //    std::cout << "num mesh-edges = " << mesh_edges.size() << std::endl;


    // group mesh-faces into patches
    std::vector<std::vector<size_t>> patches;
    {
        timing_labels.emplace_back("patches");
        ScopedTimer<> timer("patches");
        compute_patches(edges_of_mesh_face, mesh_edges, patches);
        timings.push_back(timer.toc());
    }
    std::cout << "num patches = " << patches.size() << std::endl;

    // compute map: iso-face Id --> patch Id
    //    std::vector<size_t> patch_of_face;
    //    {
    //        timing_labels.emplace_back("face-patch map");
    //        ScopedTimer<> timer("face-patch map");
    //        patch_of_face.resize(iso_faces.size());
    //        for (size_t i = 0; i < patches.size(); i++) {
    //            for (const auto& fId : patches[i]) {
    //                patch_of_face[fId] = i;
    //            }
    //        }
    //        timings.push_back(timer.toc());
    //    }

    // group non-manifold iso-edges into chains
    std::vector<std::vector<size_t>> non_manifold_edges_of_vert;
    std::vector<std::vector<size_t>> chains;
    {
        timing_labels.emplace_back("chains");
        ScopedTimer<> timer("chains");
        non_manifold_edges_of_vert.resize(mesh_verts.size());
        // get incident non-manifold edges for iso-vertices
        for (size_t i = 0; i < mesh_edges.size(); i++) {
            if (mesh_edges[i].face_edge_indices.size() >
                2) { // non-manifold edge (not a boundary edge)
                // there is only one patch incident to a boundary edge,
                // so there is no need to figure out the "order" of patches around a boundary
                // edge
                non_manifold_edges_of_vert[mesh_edges[i].v1].push_back(i);
                non_manifold_edges_of_vert[mesh_edges[i].v2].push_back(i);
            }
        }
        // group non-manifold iso-edges into chains
        compute_chains(mesh_edges, non_manifold_edges_of_vert, chains);
        timings.push_back(timer.toc());
    }
    std::cout << "num chains = " << chains.size() << std::endl;


    //
    if (!args.timing_only) {
        std::vector<std::vector<std::pair<size_t, int>>> half_patch_list;
        std::vector<std::vector<size_t>> shells;
        std::vector<std::vector<size_t>> components;
        std::vector<std::vector<size_t>> material_cells;
        save_result(output_dir + "/DC_mesh.json",
            mesh_verts,
            mesh_faces,
            patches,
            mesh_edges,
            chains,
            non_manifold_edges_of_vert,
            half_patch_list,
            shells,
            components,
            material_cells);
        save_result_msh_DC(output_dir + "/DC_mesh",
            mesh_verts,
            mesh_faces,
            patches,
            mesh_edges,
            chains,
            non_manifold_edges_of_vert);
    }
    save_timings(output_dir + "/timings.json", timing_labels, timings);

    return 0;
}