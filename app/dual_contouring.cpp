//
// Created by Charles Du on 1/19/22.
//
#include "ScopedTimer.h"
#include "implicit_arrangement_util.h"

#include <CLI/CLI.hpp>
#include <Eigen/Core>

int main(int argc, const char* argv[])
{
    struct
    {
        std::string config_file;
        bool timing_only = false;
    } args;
    CLI::App app{"Dual Contouring Command Line"};
    app.add_option("config_file", args.config_file, "Configuration file")->required();
    app.add_option("-T,--timing-only", args.timing_only, "Record timing without output result");
    CLI11_PARSE(app, argc, argv);

    // parse configure file
    std::string tet_mesh_file;
    std::string sphere_file;
    std::string output_dir;
    bool b_place_holder;
    bool use_topo_ray_shooting = true;
    std::array<double, 3> bbox_min, bbox_max;
    parse_config_file(args.config_file,
        tet_mesh_file,
        sphere_file,
        output_dir,
        b_place_holder,
        b_place_holder,
        use_topo_ray_shooting,
        b_place_holder,
        bbox_min,
        bbox_max);

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


    // load implicit function values, or evaluate
    std::vector<Sphere> spheres;
    load_spheres(sphere_file, spheres);
    size_t n_func = spheres.size();

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> funcVals;
    {
        timing_labels.emplace_back("func values");
        ScopedTimer<> timer("func values");
        funcVals.resize(n_pts, n_func);
        for (Eigen::Index i = 0; i < n_pts; i++) {
            const auto& p = pts[i];
            for (Eigen::Index j = 0; j < n_func; j++) {
                funcVals(i, j) = sphere_function(spheres[j].first, spheres[j].second, p);
            }
        }
        timings.push_back(timer.toc());
    }

    // function signs at vertices
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> funcSigns;
    {
        timing_labels.emplace_back("func signs");
        ScopedTimer<> timer("func signs");
        funcSigns.resize(n_pts, n_func);
        for (Eigen::Index i = 0; i < n_pts; i++) {
            for (Eigen::Index j = 0; j < n_func; j++) {
                funcSigns(i, j) = (funcVals(i, j) > 0);
            }
        }
        timings.push_back(timer.toc());
    }

    //    size_t num_intersecting_tet = 0;
    size_t num_1_func = 0;
    size_t num_2_func = 0;
    size_t num_more_func = 0;
    std::vector<size_t> func_in_tet;
    std::vector<size_t> start_index_of_tet;
    {
        timing_labels.emplace_back("filter");
        ScopedTimer<> timer("filter(CRS vector)");
        func_in_tet.reserve(n_tets);
        start_index_of_tet.reserve(n_tets + 1);
        start_index_of_tet.push_back(0);
        int pos_count;
        int neg_count;
        int func_count;
        for (Eigen::Index i = 0; i < n_tets; i++) {
            func_count = 0;
            for (Eigen::Index j = 0; j < n_func; j++) {
                pos_count = 0;
                neg_count = 0;
                for (size_t& vId : tets[i]) {
                    if (funcSigns(vId, j)) {
                        ++pos_count;
                    } else {
                        ++neg_count;
                    }
                }
                // tets[i].size() == 4
                if (pos_count < 4 && neg_count < 4) {
                    func_in_tet.push_back(j);
                    ++func_count;
                }
            }
            //            if (func_in_tet.size() > start_index_of_tet.back()) {
            //                ++num_intersecting_tet;
            //            }
            switch (func_count) {
            case 0: break;
            case 1: ++num_1_func; break;
            case 2: ++num_2_func; break;
            default: ++num_more_func; break;
            }
            start_index_of_tet.push_back(func_in_tet.size());
        }
        timings.push_back(timer.toc());
    }
    //    std::cout << "num_intersecting_tet = " << num_intersecting_tet << std::endl;

    // dual contouring
    std::vector<std::array<double, 3>> mesh_verts;
    std::vector<std::array<size_t, 3>> mesh_tris;
    {
        timing_labels.emplace_back("dual contouring");
        ScopedTimer<> timer("dual contouring");
        tet_dual_contouring(num_1_func,
            num_2_func,
            num_more_func,
            pts,
            tets,
            funcVals,
            funcSigns,
            func_in_tet,
            start_index_of_tet,
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
        std::vector<std::vector<size_t>> arrangement_cells;
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
            arrangement_cells);
    }
    save_timings(output_dir + "/timings.json", timing_labels, timings);


    return 0;
}