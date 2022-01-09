#include <simplicial_arrangement/lookup_table.h>
#include <simplicial_arrangement/simplicial_arrangement.h>

#include <chrono>
#include <iostream>
#include <queue>
#include "ScopedTimer.h"

#include <CLI/CLI.hpp>
#include <Eigen/Core>


#include "implicit_arrangement_util.h"

typedef std::chrono::duration<double> Time_duration;

using namespace simplicial_arrangement;


int main(int argc, const char* argv[])
{
    struct
    {
        std::string config_file;
        bool timing_only = false;
    } args;
    CLI::App app{"Implicit Arrangement Command Line"};
    app.add_option("config_file", args.config_file, "Configuration file")->required();
    app.add_option("-T,--timing-only", args.timing_only, "Record timing without output result");
    CLI11_PARSE(app, argc, argv);

    // load lookup table
    std::cout << "load table ..." << std::endl;
    bool loaded = load_lookup_table();
    if (loaded) {
        std::cout << "loading finished." << std::endl;
    } else {
        std::cout << "loading failed." << std::endl;
        return -1;
    }

    // record timings
    std::vector<std::string> timing_labels;
    std::vector<double> timings;

    // parse configure file
    std::string tet_mesh_file;
    std::string sphere_file;
    std::string output_dir;
    bool use_2func_lookup;
    parse_config_file(args.config_file, tet_mesh_file, sphere_file, output_dir, use_2func_lookup);

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
        funcVals.resize(n_func, n_pts);
        for (Eigen::Index i = 0; i < n_func; i++) {
            for (Eigen::Index j = 0; j < n_pts; j++) {
                funcVals(i, j) = sphere_function(spheres[i].first, spheres[i].second, pts[j]);
            }
        }
        timings.push_back(timer.toc());
    }


    // function signs at vertices
    Eigen::MatrixXi funcSigns;
    {
        timing_labels.emplace_back("func signs");
        ScopedTimer<> timer("func signs");
        funcSigns.resize(n_func, n_pts);
        for (Eigen::Index i = 0; i < n_func; i++) {
            for (Eigen::Index j = 0; j < n_pts; j++) {
                funcSigns(i, j) = sign(funcVals(i, j));
            }
        }
        timings.push_back(timer.toc());
    }
    //    {
    //        ScopedTimer<> timer("compute function signs at vertices");
    //        funcSigns = funcVals.unaryExpr([](double x) { return sign(x); });
    //    }

    size_t num_intersecting_tet = 0;
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
        for (Eigen::Index i = 0; i < n_tets; i++) {
            for (Eigen::Index j = 0; j < n_func; j++) {
                pos_count = 0;
                neg_count = 0;
                for (size_t& vId : tets[i]) {
                    if (funcSigns(j, vId) == 1) {
                        pos_count += 1;
                    } else if (funcSigns(j, vId) == -1) {
                        neg_count += 1;
                    }
                }
                // tets[i].size() == 4
                if (pos_count < 4 && neg_count < 4) {
                    func_in_tet.push_back(j);
                }
            }
            if (func_in_tet.size() > start_index_of_tet.back()) {
                ++num_intersecting_tet;
            }
            start_index_of_tet.push_back(func_in_tet.size());
        }
        timings.push_back(timer.toc());
    }
    std::cout << "num_intersecting_tet = " << num_intersecting_tet << std::endl;


    // compute arrangement in each tet (iterative plane cut)
    std::vector<Arrangement<3>> cut_results;
    std::vector<size_t> cut_result_index;
    //
    Time_duration time_1_func = Time_duration::zero();
    Time_duration time_2_func = Time_duration::zero();
    Time_duration time_more_func = Time_duration::zero();
    //
    {
        timing_labels.emplace_back("simp_arr(other)");
        ScopedTimer<> timer("simp_arr");
        cut_results.reserve(num_intersecting_tet);
        cut_result_index.reserve(n_tets);
        size_t start_index;
        size_t num_func;
        for (size_t i = 0; i < tets.size(); i++) {
            start_index = start_index_of_tet[i];
            num_func = start_index_of_tet[i + 1] - start_index;
            if (num_func == 0) {
                cut_result_index.push_back(Arrangement<3>::None);
                continue;
            }
            std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
            //
            size_t v1 = tets[i][0];
            size_t v2 = tets[i][1];
            size_t v3 = tets[i][2];
            size_t v4 = tets[i][3];
            std::vector<Plane<double, 3>> planes(num_func);
            for (size_t j = 0; j < num_func; j++) {
                size_t f_id = func_in_tet[start_index + j];
                planes[j] = {
                    funcVals(f_id, v1), funcVals(f_id, v2), funcVals(f_id, v3), funcVals(f_id, v4)};
            }
            //
            if (!use_2func_lookup && num_func == 2) {
                cut_result_index.push_back(cut_results.size());
                disable_lookup_table();
                cut_results.emplace_back(std::move(compute_arrangement(planes)));
                enable_lookup_table();
            } else {
                cut_result_index.push_back(cut_results.size());
                cut_results.emplace_back(std::move(compute_arrangement(planes)));
            }

            //
            std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
            switch (num_func) {
            case 1: time_1_func += std::chrono::duration_cast<Time_duration>(t2 - t1); break;
            case 2: time_2_func += std::chrono::duration_cast<Time_duration>(t2 - t1); break;
            default: time_more_func += std::chrono::duration_cast<Time_duration>(t2 - t1); break;
            }
        }
        timings.push_back(
            timer.toc() - time_1_func.count() - time_2_func.count() - time_more_func.count());
    }
    timing_labels.emplace_back("simp_arr(1 func)");
    timings.push_back(time_1_func.count());
    timing_labels.emplace_back("simp_arr(2 func)");
    timings.push_back(time_2_func.count());
    timing_labels.emplace_back("simp_arr(>=3 func)");
    timings.push_back(time_more_func.count());
    std::cout << " -- [simp_arr(1 func)]: " << time_1_func.count() << " s" << std::endl;
    std::cout << " -- [simp_arr(2 func)]: " << time_2_func.count() << " s" << std::endl;
    std::cout << " -- [simp_arr(>=3 func)]: " << time_more_func.count() << " s" << std::endl;

    // extract arrangement mesh
    std::vector<IsoVert> iso_verts;
    std::vector<IsoFace> iso_faces;
    {
        timing_labels.emplace_back("extract mesh");
        ScopedTimer<> timer("extract mesh");
        extract_iso_mesh_pure(cut_results,
            cut_result_index,
            func_in_tet,
            start_index_of_tet,
            tets,
            iso_verts,
            iso_faces);
        timings.push_back(timer.toc());
    }
    // std::cout << "num iso-vertices = " << iso_verts.size() << std::endl;
    // std::cout << "num iso-faces = " << iso_faces.size() << std::endl;

    // compute xyz of iso-vertices
    std::vector<std::array<double, 3>> iso_pts;
    {
        timing_labels.emplace_back("compute xyz");
        ScopedTimer<> timer("compute xyz");
        compute_iso_vert_xyz(iso_verts, funcVals, pts, iso_pts);
        timings.push_back(timer.toc());
    }

    //  compute iso-edges and edge-face connectivity
    std::vector<IsoEdge> iso_edges;
    {
        timing_labels.emplace_back("isoEdge-face connectivity");
        ScopedTimer<> timer("isoEdge-face connectivity");
        compute_iso_edges(iso_faces, iso_edges);
        timings.push_back(timer.toc());
    }
    // std::cout << "num iso-edges = " << iso_edges.size() << std::endl;


    // group iso-faces into patches
    std::vector<std::vector<size_t>> patches;
    {
        timing_labels.emplace_back("patches");
        ScopedTimer<> timer("patches");
        compute_patches(iso_faces, iso_edges, patches);
        timings.push_back(timer.toc());
    }
    // std::cout << "num patches = " << patches.size() << std::endl;

    // compute map: iso-face Id --> patch Id
    std::vector<size_t> patch_of_face;
    {
        timing_labels.emplace_back("face-patch map");
        ScopedTimer<> timer("face-patch map");
        patch_of_face.resize(iso_faces.size());
        for (size_t i = 0; i < patches.size(); i++) {
            for (const auto& fId : patches[i]) {
                patch_of_face[fId] = i;
            }
        }
        timings.push_back(timer.toc());
    }

    std::vector<std::vector<size_t>> arrangement_cells;
    // another approach: order patches around chains

    // group non-manifold iso-edges into chains
    std::vector<std::vector<size_t>> non_manifold_edges_of_vert;
    std::vector<std::vector<size_t>> chains;
    {
        timing_labels.emplace_back("chains");
        ScopedTimer<> timer("chains");
        non_manifold_edges_of_vert.resize(iso_pts.size());
        // get incident non-manifold edges for iso-vertices
        for (size_t i = 0; i < iso_edges.size(); i++) {
            if (iso_edges[i].face_edge_indices.size() >
                2) { // non-manifold edge (not a boundary edge)
                // there is only one patch incident to a boundary edge,
                // so there is no need to figure out the "order" of patches around a boundary
                // edge
                non_manifold_edges_of_vert[iso_edges[i].v1].push_back(i);
                non_manifold_edges_of_vert[iso_edges[i].v2].push_back(i);
            }
        }
        // group non-manifold iso-edges into chains
        compute_chains(iso_edges, non_manifold_edges_of_vert, chains);
        timings.push_back(timer.toc());
    }
    // std::cout << "num chains = " << chains.size() << std::endl;


    {
        timing_labels.emplace_back("vert-tet connectivity");
        ScopedTimer<> timer("vert-tet connectivity(compact vec)");
        // compute number of tets incident to each vertex
        std::vector<int> num_incident_tets(n_pts, 0);
        for (const auto& tet : tets) {
            num_incident_tets[tet[0]] += 1;
            num_incident_tets[tet[1]] += 1;
            num_incident_tets[tet[2]] += 1;
            num_incident_tets[tet[3]] += 1;
        }
        // get index of pts in
        std::vector<int> pts_index(n_pts);
        int curr_index = 0;
        for (int i = 0; i < n_pts; ++i) {
            pts_index[i] = curr_index;
            curr_index += num_incident_tets[i];
        }
        //
        std::vector<int> size_incident_tets(n_pts, 0);
        std::vector<int> incident_tets(4 * n_tets);
        for (int i = 0; i < n_tets; ++i) {
            const auto& tet = tets[i];
            incident_tets[pts_index[tet[0]] + size_incident_tets[tet[0]]] = i;
            incident_tets[pts_index[tet[1]] + size_incident_tets[tet[1]]] = i;
            incident_tets[pts_index[tet[2]] + size_incident_tets[tet[2]]] = i;
            incident_tets[pts_index[tet[3]] + size_incident_tets[tet[3]]] = i;
            size_incident_tets[tet[0]] += 1;
            size_incident_tets[tet[1]] += 1;
            size_incident_tets[tet[2]] += 1;
            size_incident_tets[tet[3]] += 1;
        }
        timings.push_back(timer.toc());
        //            std::cout << "num_incident_tets[0] = " << num_incident_tets[0] << std::endl;
        //            std::cout << "pts_index[0] = " << pts_index[0] << std::endl;
        //            std::cout << "incident_tets[0] = " << incident_tets[0] << std::endl;
    }

    //        {
    // this approach is slower than (compact vec)
    //            timing_labels.emplace_back("vert-tet connectivity(Eigen)");
    //            ScopedTimer<> timer("vert-tet connectivity(Eigen)");
    //            Eigen::VectorXi num_incident_tets;
    //            num_incident_tets.setZero(n_pts);
    //            // compute number of tets incident to each vertex
    //            for (const auto & tet : tets) {
    //                num_incident_tets[tet[0]] += 1;
    //                num_incident_tets[tet[1]] += 1;
    //                num_incident_tets[tet[2]] += 1;
    //                num_incident_tets[tet[3]] += 1;
    //            }
    //            // get max
    //            int max_num_tets = num_incident_tets.maxCoeff();
    ////            std::cout << "max_num_tets = " << max_num_tets << std::endl;
    //            // compute indices of tets incident to each vertex
    //            Eigen::MatrixXi incident_tets_of_vert(n_pts, max_num_tets);
    //            Eigen::VectorXi size_incident_tets;
    //            size_incident_tets.setZero(n_pts);
    //            for (size_t i = 0; i < tets.size(); i++) {
    //                const auto& tet = tets[i];
    //                incident_tets_of_vert(tet[0],size_incident_tets[tet[0]]) = i;
    //                incident_tets_of_vert(tet[1],size_incident_tets[tet[1]]) = i;
    //                incident_tets_of_vert(tet[2],size_incident_tets[tet[2]]) = i;
    //                incident_tets_of_vert(tet[3],size_incident_tets[tet[3]]) = i;
    //                size_incident_tets[tet[0]] += 1;
    //                size_incident_tets[tet[1]] += 1;
    //                size_incident_tets[tet[2]] += 1;
    //                size_incident_tets[tet[3]] += 1;
    //            }
    //            timings.push_back(timer.toc());
    ////            std::cout << "incident_tets_of_vertA[x][y] = " << incident_tets_of_vert(0,0) << std::endl;
    //        }

    //        {
    // this approach cost more memory than (compact vec)
    //            timing_labels.emplace_back("vert-tet connectivity(flat 2D vector)");
    //            ScopedTimer<> timer("vert-tet connectivity(flat 2D vector)");
    //            std::vector<int> num_incident_tetsA(n_pts,0);
    //            // compute number of tets incident to each vertex
    //            for (const auto & tet : tets) {
    //                num_incident_tetsA[tet[0]] += 1;
    //                num_incident_tetsA[tet[1]] += 1;
    //                num_incident_tetsA[tet[2]] += 1;
    //                num_incident_tetsA[tet[3]] += 1;
    //            }
    //            // get max
    //            int max_num_tets = *(std::max_element(num_incident_tetsA.begin(),
    //            num_incident_tetsA.end()));
    ////            std::cout << "max_num_tets = " << max_num_tets << std::endl;
    //            // compute indices of tets incident to each vertex
    //            std::vector<int> incident_tets_of_vertA(n_pts * max_num_tets);
    //            std::vector<int> size_incident_tetsA(n_pts,0);
    //            for (size_t i = 0; i < tets.size(); i++) {
    //                const auto& tet = tets[i];
    //                incident_tets_of_vertA[tet[0]*max_num_tets + size_incident_tetsA[tet[0]]] = i;
    //                incident_tets_of_vertA[tet[1]*max_num_tets + size_incident_tetsA[tet[1]]] = i;
    //                incident_tets_of_vertA[tet[2]*max_num_tets + size_incident_tetsA[tet[2]]] = i;
    //                incident_tets_of_vertA[tet[3]*max_num_tets + size_incident_tetsA[tet[3]]] = i;
    //                size_incident_tetsA[tet[0]] += 1;
    //                size_incident_tetsA[tet[1]] += 1;
    //                size_incident_tetsA[tet[2]] += 1;
    //                size_incident_tetsA[tet[3]] += 1;
    //            }
    ////            std::cout << "incident_tets_of_vertA[0] = " << incident_tets_of_vertA[0] << std::endl;
    //            timings.push_back(timer.toc());
    //        }


    // compute order of patches around chains
    // pair<size_t, int> : pair (iso-face index, iso-face orientation)
    std::vector<std::vector<std::pair<size_t, int>>> half_faces_list;
    std::vector<std::vector<std::pair<size_t, int>>> half_patch_list;
    {
        timing_labels.emplace_back("order patches around chains");
        ScopedTimer<> timer("order patches around chains");
        half_faces_list.resize(chains.size());
        half_patch_list.resize(half_faces_list.size());
        // pick representative iso-edge from each chain
        std::vector<size_t> chain_representatives(chains.size());
        for (size_t i = 0; i < chains.size(); i++) {
            chain_representatives[i] = chains[i][0];
        }
        // order iso-faces incident to each representative iso-edge
        for (size_t i = 0; i < chain_representatives.size(); i++) {
            // first try: assume each representative edge is in the interior of a tetrahedron
            const auto& iso_edge = iso_edges[chain_representatives[i]];
            auto iso_face_id = iso_edge.face_edge_indices[0].first;
            auto tet_id = iso_faces[iso_face_id].tet_face_indices[0].first;
            compute_face_order_in_one_tet(
                cut_results[cut_result_index[tet_id]], iso_faces, iso_edge, half_faces_list[i]);
        }
        // replace iso-face indices by patch indices
        for (size_t i = 0; i < half_faces_list.size(); i++) {
            half_patch_list[i].resize(half_faces_list[i].size());
            for (size_t j = 0; j < half_faces_list[i].size(); j++) {
                half_patch_list[i][j] = std::make_pair(
                    patch_of_face[half_faces_list[i][j].first], half_faces_list[i][j].second);
            }
        }
        timings.push_back(timer.toc());
    }

    // group patches into arrangement cells
    // each cell is a represented as a list of patch indices
    {
        timing_labels.emplace_back("arrangement cells");
        ScopedTimer<> timer("arrangement cells");
        compute_arrangement_cells(patches.size(), half_patch_list, arrangement_cells);
        timings.push_back(timer.toc());
    }
    std::cout << "num arrangement cells = " << arrangement_cells.size() << std::endl;

    // test: export iso-mesh, patches, chains
    if (!args.timing_only) {
        save_result(output_dir + "/iso_mesh.json",
            iso_pts,
            iso_faces,
            patches,
            iso_edges,
            chains,
            non_manifold_edges_of_vert,
            half_patch_list,
            arrangement_cells);
    }
    // test: export timings
    save_timings(output_dir + "/timings.json", timing_labels, timings);

    return 0;
}
