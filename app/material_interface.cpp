//
// Created by Charles Du on 1/15/22.
//

#include <simplicial_arrangement/lookup_table.h>
//#include <simplicial_arrangement/simplicial_arrangement.h>
#include <simplicial_arrangement/material_interface.h>

#include <chrono>
#include <iostream>
#include <queue>
#include "ScopedTimer.h"

#include <absl/container/flat_hash_map.h>
#include <absl/container/flat_hash_set.h>
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
        bool robust_test = false;
    } args;
    CLI::App app{"Material Interface Command Line"};
    app.add_option("config_file", args.config_file, "Configuration file")->required();
    app.add_option("-T,--timing-only", args.timing_only, "Record timing without output result");
    app.add_option("-R,--robust-test",args.robust_test, "Perform robustness test");
    CLI11_PARSE(app, argc, argv);

    // parse configure file
    std::string tet_mesh_file;
    std::string material_file;
    std::string output_dir;
    bool use_lookup = true;
    bool use_3func_lookup = true;
    bool use_topo_ray_shooting = true;
    parse_config_file_MI(args.config_file, tet_mesh_file, material_file, output_dir,
        use_lookup, use_3func_lookup, use_topo_ray_shooting);
    //    std::string config_path = args.config_file.substr(0, args.config_file.find_last_of('/'));
    //    std::cout << "config path: " << config_path << std::endl;
    //    tet_mesh_file = config_path + "/" + tet_mesh_file;
    //    material_file = config_path + "/" + material_file;
    if (use_lookup) {
        // load lookup table
        std::cout << "load table ..." << std::endl;
        bool loaded = load_lookup_table(simplicial_arrangement::MATERIAL_INTERFACE);
        if (loaded) {
            std::cout << "loading finished." << std::endl;
        } else {
            std::cout << "loading failed." << std::endl;
            return -1;
        }
    } else {
        use_3func_lookup = false;
    }

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
//    std::cout << "pts[100] = (" << pts[100][0] << "," << pts[100][1] << "," << pts[100][2] << ")" << std::endl;
//    std::cout << "tets[100] = (" << tets[100][0] << "," << tets[100][1] << "," << tets[100][2] << "," <<
//        tets[100][3] << ")" << std::endl;

    // robustness test
    if (args.robust_test) {
        // load robustness test spheres
        std::vector<Sphere> spheres;
        load_spheres(material_file, spheres);
        size_t n_shift_sphere = 4;
        size_t n_test = spheres.size() / n_shift_sphere;
        // the last function is 0
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> funcVals;
        size_t n_func = n_shift_sphere + 1;
        funcVals.setZero(n_pts, n_func);
        //
        std::vector<size_t> highest_material_at_vert;
        highest_material_at_vert.reserve(n_pts);
        std::vector<bool> is_degenerate_vertex;
        is_degenerate_vertex.reserve(n_pts);
        bool found_degenerate_vertex = false;
        absl::flat_hash_map<size_t, std::vector<size_t>> highest_materials_at_vert;
        //
        size_t num_intersecting_tet = 0;
        std::vector<size_t> material_in_tet;
        std::vector<size_t> start_index_of_tet;
        material_in_tet.reserve(n_tets);
        start_index_of_tet.reserve(n_tets + 1);
        //
        MaterialInterface<3> cut_result, cut_result2;
        std::vector<Material<double, 3>> materials;
        std::vector<Material<double, 3>> materials_reverse; // materials in reverse order
        materials.reserve(3);
        materials_reverse.reserve(3);
        //
        size_t type3_count = 0; // normal order run, reverse order crash
        size_t type2_count = 0; // normal order crash
        size_t type1_count = 0; // normal and reverse order don't crash, but are different
        //
        for (size_t iter = 0; iter < n_test; ++iter) {
            std::cout << "-------- test " << iter << " -----------" << std::endl;
            // load implicit functions
            for (Eigen::Index i = 0; i < n_pts; ++i) {
                const auto& p = pts[i];
                for (Eigen::Index j = 0; j < n_shift_sphere; j++) {
                    funcVals(i, j) = compute_sphere_distance(
                        spheres[iter * n_shift_sphere + j].first,
                        spheres[iter * n_shift_sphere + j].second, p);
                }
            }
            // highest material at vertices
            {
                timing_labels.emplace_back("highest func");
                ScopedTimer<> timer("highest func");
                is_degenerate_vertex.clear();
                is_degenerate_vertex.resize(n_pts, false);
                highest_material_at_vert.clear();
                highest_materials_at_vert.clear();
                found_degenerate_vertex = false;
                for (Eigen::Index i = 0; i < n_pts; i++) {
                    double max = funcVals(i, 0);
                    size_t max_id = 0;
                    size_t max_count = 1;
                    for (Eigen::Index j = 1; j < n_func; j++) {
                        if (funcVals(i, j) > max) {
                            max = funcVals(i, j);
                            max_id = j;
                            max_count = 1;
                        } else if (funcVals(i, j) == max) {
                            ++max_count;
                        }
                    }
                    highest_material_at_vert.push_back(max_id);
                    //
                    if (max_count > 1) {
                        is_degenerate_vertex[i] = true;
                        found_degenerate_vertex = true;
                        auto& materials = highest_materials_at_vert[i];
                        materials.reserve(max_count);
                        for (Eigen::Index j = 0; j < n_func; j++) {
                            if (funcVals(i, j) == max) {
                                materials.push_back(j);
                            }
                        }
                    }
                }
                timings.push_back(timer.toc());
            }
            // filter
            {
                timing_labels.emplace_back("filter");
                ScopedTimer<> timer("filter(CRS vector)");
                material_in_tet.clear();
                start_index_of_tet.clear();
                start_index_of_tet.push_back(0);
                std::set<size_t> materials;
                std::array<double, 4> min_h;
                for (size_t i = 0; i < n_tets; ++i) {
                    const auto& tet = tets[i];
                    // find high materials
                    materials.clear();
                    for (size_t j = 0; j < 4; ++j) {
                        if (is_degenerate_vertex[tet[j]]) {
                            const auto& ms = highest_materials_at_vert[tet[j]];
                            materials.insert(ms.begin(), ms.end());
                        } else {
                            materials.insert(highest_material_at_vert[tet[j]]);
                        }
                    }
                    // if only one high material, there is no material interface
                    if (materials.size() < 2) { // no material interface
                        start_index_of_tet.push_back(material_in_tet.size());
                        continue;
                    }
                    // find min of high materials
                    min_h.fill(std::numeric_limits<double>::max());
                    for (auto it = materials.begin(); it != materials.end(); it++) {
                        for (size_t j = 0; j < 4; ++j) {
                            if (funcVals(tet[j], *it) < min_h[j]) {
                                min_h[j] = funcVals(tet[j], *it);
                            }
                        }
                    }
                    // find materials greater than at least two mins of high materials
                    size_t greater_count;
                    for (size_t j = 0; j < n_func; ++j) {
                        greater_count = 0;
                        for (size_t k = 0; k < 4; ++k) {
                            if (funcVals(tet[k], j) > min_h[k]) {
                                ++greater_count;
                            }
                        }
                        if (greater_count > 1) {
                            materials.insert(j);
                        }
                    }
                    //
                    ++num_intersecting_tet;
                    material_in_tet.insert(
                        material_in_tet.end(), materials.begin(), materials.end());
                    start_index_of_tet.push_back(material_in_tet.size());
                }
                timings.push_back(timer.toc());
            }
            std::cout << "num_intersecting_tet = " << num_intersecting_tet << std::endl;
            // compute material interface in each tet
            //            std::vector<MaterialInterface<3>> cut_results;
            //            std::vector<size_t> cut_result_index;
            //
            Time_duration time_2_func = Time_duration::zero();
            Time_duration time_3_func = Time_duration::zero();
            Time_duration time_more_func = Time_duration::zero();
            size_t num_2_func = 0;
            size_t num_3_func = 0;
            size_t num_more_func = 0;
            //
            bool is_type3 = false;
            bool is_type2 = false;
            bool is_type1 = false;
            //
            {
                timing_labels.emplace_back("MI(other)");
                ScopedTimer<> timer("material interface in tets");
                if (args.robust_test) {
                    size_t start_index;
                    size_t num_func;
                    for (size_t i = 0; i < tets.size(); i++) {
                        start_index = start_index_of_tet[i];
                        num_func = start_index_of_tet[i + 1] - start_index;
                        if (num_func == 0) {
                            continue;
                        }
                        std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
                        //
                        size_t v1 = tets[i][0];
                        size_t v2 = tets[i][1];
                        size_t v3 = tets[i][2];
                        size_t v4 = tets[i][3];
                        materials.clear();
                        for (size_t j = 0; j < num_func; j++) {
                            size_t f_id = material_in_tet[start_index + j];
                            materials.emplace_back();
                            auto& material = materials.back();
                            material[0] = funcVals(v1, f_id);
                            material[1] = funcVals(v2, f_id);
                            material[2] = funcVals(v3, f_id);
                            material[3] = funcVals(v4, f_id);
                        }
                        // reverse material order
                        materials_reverse.clear();
                        for (size_t j = 0; j < num_func; ++j) {
                            materials_reverse.emplace_back(materials[num_func - 1 - j]);
                        }
                        //
                        bool crashed = false;
                        if (!use_3func_lookup && num_func == 3) {
                            disable_lookup_table();
                            try {
                                cut_result = compute_material_interface(materials);
                            } catch (std::runtime_error& e) {
                                crashed = true;
                                is_type2 = true;
                                break;
                            }
                            if (!crashed) {
                                try {
                                    cut_result2 = compute_material_interface(materials_reverse);
                                } catch (std::runtime_error& e) {
                                    crashed = true;
//                                    is_type2 = true;
                                    is_type3 = true;
//                                    break;
                                }
                                if (!crashed) {
                                    if (cut_result.vertices.size() !=
                                            cut_result2.vertices.size() ||
                                        cut_result.faces.size() !=
                                            cut_result2.faces.size() ||
                                        cut_result.cells.size() !=
                                            cut_result2.cells.size()) {
                                        // inconsistent results
                                        is_type1 = true;
//                                        break;
                                    }
                                }
                            }
                            enable_lookup_table();
                        } else {
                            try {
                                cut_result = compute_material_interface(materials);
                            } catch (std::runtime_error& e) {
                                crashed = true;
                                is_type2 = true;
                                break;
                            }
                            if (!crashed) {
                                try {
                                    cut_result2 = compute_material_interface(materials_reverse);
                                } catch (std::runtime_error& e) {
                                    crashed = true;
//                                    is_type2 = true;
                                    is_type3 = true;
//                                    break;
                                }
                                if (!crashed) {
                                    if (cut_result.vertices.size() !=
                                            cut_result2.vertices.size() ||
                                        cut_result.faces.size() !=
                                            cut_result2.faces.size() ||
                                        cut_result.cells.size() !=
                                            cut_result2.cells.size()) {
                                        // inconsistent results
                                        is_type1 = true;
                                    }
                                }
                            }
                        }
                        //
                        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
                        switch (num_func) {
                        case 2:
                            time_2_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                            ++num_2_func;
                            break;
                        case 3:
                            time_3_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                            ++num_3_func;
                            break;
                        default:
                            time_more_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                            ++num_more_func;
                            break;
                        }
                    }
                }
                timings.push_back(timer.toc() - time_2_func.count() - time_3_func.count() -
                                  time_more_func.count());
            }
            timing_labels.emplace_back("MI(2 func)");
            timings.push_back(time_2_func.count());
            timing_labels.emplace_back("MI(3 func)");
            timings.push_back(time_3_func.count());
            timing_labels.emplace_back("MI(>=4 func)");
            timings.push_back(time_more_func.count());
            std::cout << " -- [MI(2 func)]: " << time_2_func.count() << " s" << std::endl;
            std::cout << " -- [MI(3 func)]: " << time_3_func.count() << " s" << std::endl;
            std::cout << " -- [MI(>=4 func)]: " << time_more_func.count() << " s" << std::endl;

            if (is_type2) {
                std::cout << "type 2 failure." << std::endl;
                ++type2_count;
            } else if (is_type3) {
                std::cout << "type 3 failure." << std::endl;
                ++type3_count;
            } else if (is_type1) {
                std::cout << "type 1 failure." << std::endl;
                ++type1_count;
            }

            std::cout << "=======================" << std::endl;
            std::cout << "total: " << iter + 1 << std::endl;
            std::cout << "type 1: " << type1_count << std::endl;
            std::cout << "type 2: " << type2_count << std::endl;
            std::cout << "type 3: " << type3_count << std::endl;
        }
        return 0;

    }

    // load implicit functions and compute function values at vertices
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> funcVals;
    if (load_functions(material_file, pts, funcVals)) {
        std::cout << "function loading finished." << std::endl;
    } else {
        std::cout << "function loading failed." << std::endl;
        return -2;
    }
    size_t n_func = funcVals.cols();
    std::cout << "n_func = " << n_func << std::endl;


    // highest material at vertices
    std::vector<size_t> highest_material_at_vert;
    // degenerate vertex: more than one highest material, i.e. material interface passes the vertex
    std::vector<bool> is_degenerate_vertex;
    bool found_degenerate_vertex = false;
    absl::flat_hash_map<size_t, std::vector<size_t>> highest_materials_at_vert;
    {
        timing_labels.emplace_back("highest func");
        ScopedTimer<> timer("highest func");
        is_degenerate_vertex.resize(n_pts, false);
        highest_material_at_vert.reserve(n_pts);
        for (Eigen::Index i = 0; i < n_pts; i++) {
            double max = funcVals(i, 0);
            size_t max_id = 0;
            size_t max_count = 1;
            for (Eigen::Index j = 1; j < n_func; j++) {
                if (funcVals(i,j) > max) {
                    max = funcVals(i,j);
                    max_id = j;
                    max_count = 1;
                } else if (funcVals(i,j) == max) {
                    ++max_count;
                }
            }
            highest_material_at_vert.push_back(max_id);
            //
            if (max_count > 1) {
                is_degenerate_vertex[i] = true;
                found_degenerate_vertex = true;
                auto& materials = highest_materials_at_vert[i];
                materials.reserve(max_count);
                for (Eigen::Index j = 0; j < n_func; j++) {
                    if (funcVals(i,j) == max) {
                        materials.push_back(j);
                    }
                }
            }
        }
        timings.push_back(timer.toc());
    }

    // filter relevant materials in each tet
    // a tet is non-empty if there are some material interface in it
    size_t num_intersecting_tet = 0;
    std::vector<size_t> material_in_tet;
    std::vector<size_t> start_index_of_tet;
    {
        timing_labels.emplace_back("filter");
        ScopedTimer<> timer("filter(CRS vector)");
        material_in_tet.reserve(n_tets);
        start_index_of_tet.reserve(n_tets + 1);
        start_index_of_tet.push_back(0);
        std::set<size_t> materials;
        std::array<double, 4> min_h;
        for (size_t i = 0; i < n_tets; ++i) {
            const auto& tet = tets[i];
            // find high materials
            materials.clear();
            for (size_t j = 0; j < 4; ++j) {
                if (is_degenerate_vertex[tet[j]]) {
                    const auto& ms = highest_materials_at_vert[tet[j]];
                    materials.insert(ms.begin(), ms.end());
                } else {
                    materials.insert(highest_material_at_vert[tet[j]]);
                }
            }
            // if only one high material, there is no material interface
            if (materials.size() < 2) {  // no material interface
                start_index_of_tet.push_back(material_in_tet.size());
                continue;
            }
            // find min of high materials
            min_h.fill(std::numeric_limits<double>::max());
            for(auto it = materials.begin(); it != materials.end(); it++) {
                for (size_t j = 0; j < 4; ++j) {
                    if (funcVals(tet[j], *it) < min_h[j]) {
                        min_h[j] = funcVals(tet[j], *it);
                    }
                }
            }
            // find materials greater than at least two mins of high materials
            size_t greater_count;
            for (size_t j = 0; j < n_func; ++j) {
                greater_count = 0;
                for (size_t k = 0; k < 4; ++k) {
                    if (funcVals(tet[k],j) > min_h[k]) {
                        ++greater_count;
                    }
                }
                if (greater_count > 1) {
                    materials.insert(j);
                }
            }
            //
            ++num_intersecting_tet;
            material_in_tet.insert(material_in_tet.end(), materials.begin(), materials.end());
            start_index_of_tet.push_back(material_in_tet.size());
        }
        timings.push_back(timer.toc());
    }
    std::cout << "num_intersecting_tet = " << num_intersecting_tet << std::endl;


    // compute material interface in each tet
    std::vector<MaterialInterface<3>> cut_results;
    std::vector<size_t> cut_result_index;
    //
    Time_duration time_2_func = Time_duration::zero();
    Time_duration time_3_func = Time_duration::zero();
    Time_duration time_more_func = Time_duration::zero();
    size_t num_2_func = 0;
    size_t num_3_func = 0;
    size_t num_more_func = 0;
    //
    size_t type2_count = 0;
    size_t type1_count = 0;
    //
    {
        timing_labels.emplace_back("MI(other)");
        ScopedTimer<> timer("material interface in tets");
        if (args.robust_test) {
            cut_results.reserve(num_intersecting_tet);
            cut_result_index.reserve(n_tets);
            std::vector<MaterialInterface<3>> cut_results2; // cut_results when reverse function order
            size_t start_index;
            size_t num_func;
            std::vector<Material<double, 3>> materials;
            std::vector<Material<double, 3>> materials_reverse; // materials in reverse order
            materials.reserve(3);
            materials_reverse.reserve(3);
            for (size_t i = 0; i < tets.size(); i++) {
                start_index = start_index_of_tet[i];
                num_func = start_index_of_tet[i + 1] - start_index;
                if (num_func == 0) {
                    cut_result_index.push_back(MaterialInterface<3>::None);
                    continue;
                }
                std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
                //
                size_t v1 = tets[i][0];
                size_t v2 = tets[i][1];
                size_t v3 = tets[i][2];
                size_t v4 = tets[i][3];
                materials.clear();
                for (size_t j = 0; j < num_func; j++) {
                    size_t f_id = material_in_tet[start_index + j];
                    materials.emplace_back();
                    auto& material = materials.back();
                    material[0] = funcVals(v1, f_id);
                    material[1] = funcVals(v2, f_id);
                    material[2] = funcVals(v3, f_id);
                    material[3] = funcVals(v4, f_id);
                }
                // reverse material order
                materials_reverse.clear();
                for (size_t j = 0; j < num_func; ++j) {
                    materials_reverse.emplace_back(materials[num_func -1 -j]);
                }
                //
                bool crashed = false;
                if (!use_3func_lookup && num_func == 3) {
                    cut_result_index.push_back(cut_results.size());
                    disable_lookup_table();
                    try {
                        cut_results.emplace_back(std::move(compute_material_interface(materials)));
                    } catch (std::runtime_error& e) {
                        crashed = true;
                        type2_count++;
                    }
                    if (!crashed) {
                        try {
                            cut_results2.emplace_back(std::move(compute_material_interface(materials_reverse)));
                        } catch (std::runtime_error& e) {
                            crashed = true;
                            type2_count++;
                        }
                        if (!crashed) {
                            const auto& last_cut_result = cut_results.back();
                            const auto& last_cut_result2 = cut_results2.back();
                            if (last_cut_result.vertices.size() != last_cut_result2.vertices.size() ||
                                last_cut_result.faces.size() != last_cut_result2.faces.size() ||
                                last_cut_result.cells.size() != last_cut_result2.cells.size()) {
                                // inconsistent results
                                type1_count++;
                            }
                        }
                    }
                    enable_lookup_table();
                } else {
                    cut_result_index.push_back(cut_results.size());
                    try {
                        cut_results.emplace_back(std::move(compute_material_interface(materials)));
                    } catch (std::runtime_error& e) {
                        crashed = true;
                        type2_count++;
                    }
                    if (!crashed) {
                        try {
                            cut_results2.emplace_back(std::move(compute_material_interface(materials_reverse)));
                        } catch (std::runtime_error& e) {
                            crashed = true;
                            type2_count++;
                        }
                        if (!crashed) {
                            const auto& last_cut_result = cut_results.back();
                            const auto& last_cut_result2 = cut_results2.back();
                            if (last_cut_result.vertices.size() != last_cut_result2.vertices.size() ||
                                last_cut_result.faces.size() != last_cut_result2.faces.size() ||
                                last_cut_result.cells.size() != last_cut_result2.cells.size()) {
                                // inconsistent results
                                type1_count++;
                            }
                        }
                    }
                }
                //
                std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
                switch (num_func) {
                case 2:
                    time_2_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                    ++num_2_func;
                    break;
                case 3:
                    time_3_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                    ++num_3_func;
                    break;
                default:
                    time_more_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                    ++num_more_func;
                    break;
                }
            }
        } else {  // not performing robustness test
            cut_results.reserve(num_intersecting_tet);
            cut_result_index.reserve(n_tets);
            size_t start_index;
            size_t num_func;
            std::vector<Material<double, 3>> materials;
            materials.reserve(3);
            try {
                for (size_t i = 0; i < tets.size(); i++) {
                    start_index = start_index_of_tet[i];
                    num_func = start_index_of_tet[i + 1] - start_index;
                    if (num_func == 0) {
                        cut_result_index.push_back(MaterialInterface<3>::None);
                        continue;
                    }
                    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
                    //
                    size_t v1 = tets[i][0];
                    size_t v2 = tets[i][1];
                    size_t v3 = tets[i][2];
                    size_t v4 = tets[i][3];
                    materials.clear();
                    for (size_t j = 0; j < num_func; j++) {
                        size_t f_id = material_in_tet[start_index + j];
                        materials.emplace_back();
                        auto& material = materials.back();
                        material[0] = funcVals(v1, f_id);
                        material[1] = funcVals(v2, f_id);
                        material[2] = funcVals(v3, f_id);
                        material[3] = funcVals(v4, f_id);
                    }
                    //
                    if (!use_3func_lookup && num_func == 3) {
                        cut_result_index.push_back(cut_results.size());
                        disable_lookup_table();
                        cut_results.emplace_back(std::move(compute_material_interface(materials)));
                        enable_lookup_table();
                    } else {
                        cut_result_index.push_back(cut_results.size());
                        cut_results.emplace_back(std::move(compute_material_interface(materials)));
                    }

                    //
                    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
                    switch (num_func) {
                    case 2:
                        time_2_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                        ++num_2_func;
                        break;
                    case 3:
                        time_3_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                        ++num_3_func;
                        break;
                    default:
                        time_more_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                        ++num_more_func;
                        break;
                    }
                }
            } catch (std::runtime_error& e) {
                std::cout << e.what() << std::endl;
                return -1;
            }
        }
        timings.push_back(
            timer.toc() - time_2_func.count() - time_3_func.count() - time_more_func.count());
    }
    timing_labels.emplace_back("MI(2 func)");
    timings.push_back(time_2_func.count());
    timing_labels.emplace_back("MI(3 func)");
    timings.push_back(time_3_func.count());
    timing_labels.emplace_back("MI(>=4 func)");
    timings.push_back(time_more_func.count());
    std::cout << " -- [MI(2 func)]: " << time_2_func.count() << " s" << std::endl;
    std::cout << " -- [MI(3 func)]: " << time_3_func.count() << " s" << std::endl;
    std::cout << " -- [MI(>=4 func)]: " << time_more_func.count() << " s" << std::endl;

    if (args.robust_test) {
        std::cout << "total: " << num_intersecting_tet << std::endl;
        std::cout << "type 1: " << type1_count << std::endl;
        std::cout << "type 2: " << type2_count << std::endl;
        return 0;
    }

    // extract material interface mesh
    std::vector<MI_Vert> MI_verts;
    std::vector<PolygonFace> MI_faces;
    // the following data are only needed when we use
    // the baseline nesting algorithm
    // (group simplicial cells into material cells)
    std::vector<long long> global_vId_of_tet_vert;
    std::vector<size_t> global_vId_start_index_of_tet;
    std::vector<size_t> MI_fId_of_tet_face;
    std::vector<size_t> MI_fId_start_index_of_tet;
    {
        timing_labels.emplace_back("extract mesh");
        ScopedTimer<> timer("extract mesh");
        if (use_topo_ray_shooting) {
            extract_MI_mesh_pure(num_2_func,
                num_3_func,
                num_more_func,
                cut_results,
                cut_result_index,
                material_in_tet,
                start_index_of_tet,
                tets,
                MI_verts,
                MI_faces);
        } else { // nesting algorithm: group simplicial cells into  material cells
            extract_MI_mesh(num_2_func,
                num_3_func,
                num_more_func,
                cut_results,
                cut_result_index,
                material_in_tet,
                start_index_of_tet,
                tets,
                MI_verts,
                MI_faces,
                global_vId_of_tet_vert,
                global_vId_start_index_of_tet,
                MI_fId_of_tet_face,
                MI_fId_start_index_of_tet);
        }
        timings.push_back(timer.toc());
    }
    std::cout << "num MI-vertices = " << MI_verts.size() << std::endl;
    std::cout << "num MI-faces = " << MI_faces.size() << std::endl;

    // compute xyz of MI-vertices
    std::vector<std::array<double, 3>> MI_pts;
    {
        timing_labels.emplace_back("compute xyz");
        ScopedTimer<> timer("compute xyz");
        compute_MI_vert_xyz(MI_verts, funcVals, pts, MI_pts);
        timings.push_back(timer.toc());
    }

    //  compute iso-edges and edge-face connectivity
    std::vector<Edge> MI_edges;
    std::vector<std::vector<size_t>> edges_of_MI_face;
    {
        timing_labels.emplace_back("edge-face connectivity");
        ScopedTimer<> timer("edge-face connectivity");
        compute_mesh_edges(MI_faces, edges_of_MI_face, MI_edges);
        timings.push_back(timer.toc());
    }
    // std::cout << "num iso-edges = " << iso_edges.size() << std::endl;

    // group iso-faces into patches
    std::vector<std::vector<size_t>> patches;
    {
        timing_labels.emplace_back("patches");
        ScopedTimer<> timer("patches");
        compute_patches(edges_of_MI_face, MI_edges, patches);
        timings.push_back(timer.toc());
    }
    std::cout << "num patches = " << patches.size() << std::endl;

    // compute map: MI-face Id --> patch Id
    std::vector<size_t> patch_of_face;
    {
        timing_labels.emplace_back("face-patch map");
        ScopedTimer<> timer("face-patch map");
        patch_of_face.resize(MI_faces.size());
        for (size_t i = 0; i < patches.size(); i++) {
            for (const auto& fId : patches[i]) {
                patch_of_face[fId] = i;
            }
        }
        timings.push_back(timer.toc());
    }


    // group non-manifold iso-edges into chains
    std::vector<std::vector<size_t>> non_manifold_edges_of_vert;
    std::vector<std::vector<size_t>> chains;
    {
        timing_labels.emplace_back("chains");
        ScopedTimer<> timer("chains");
        non_manifold_edges_of_vert.resize(MI_pts.size());
        // get incident non-manifold edges for MI-vertices
        for (size_t i = 0; i < MI_edges.size(); i++) {
            if (MI_edges[i].face_edge_indices.size() >
                2) { // non-manifold edge (not a boundary edge)
                // there is only one patch incident to a boundary edge,
                // so there is no need to figure out the "order" of patches around a boundary
                // edge
                non_manifold_edges_of_vert[MI_edges[i].v1].push_back(i);
                non_manifold_edges_of_vert[MI_edges[i].v2].push_back(i);
            }
        }
        // group non-manifold iso-edges into chains
        compute_chains(MI_edges, non_manifold_edges_of_vert, chains);
        timings.push_back(timer.toc());
    }
    std::cout << "num chains = " << chains.size() << std::endl;

    // vert-tet connectivity
    absl::flat_hash_map<size_t, std::vector<size_t>> incident_tets;
    {
        timing_labels.emplace_back("vert-tet connectivity");
        ScopedTimer<> timer("vert-tet connectivity(degenerate vert only)");
        if (found_degenerate_vertex) {
            for (size_t i = 0; i < n_tets; ++i) {
                const auto& tet = tets[i];
                for (size_t j = 0; j < 4; ++j) {
                    if (is_degenerate_vertex[tet[j]]) {
                        incident_tets[tet[j]].push_back(i);
                    }
                }
            }
        }
        timings.push_back(timer.toc());
    }
    //    std::cout << "found_degenerate_vertex = " << found_degenerate_vertex << std::endl;
    std::cout << "incident_tets.size() = " << incident_tets.size() << std::endl;


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
            const auto& MI_edge = MI_edges[chain_representatives[i]];
            auto MI_face_id = MI_edge.face_edge_indices[0].first;
            auto tet_id = MI_faces[MI_face_id].tet_face_indices[0].first;
            compute_face_order_in_one_tet_MI(
                cut_results[cut_result_index[tet_id]], MI_faces, MI_edge, half_faces_list[i]);
        }
        // replace MI-face indices by patch indices
        for (size_t i = 0; i < half_faces_list.size(); i++) {
            half_patch_list[i].resize(half_faces_list[i].size());
            for (size_t j = 0; j < half_faces_list[i].size(); j++) {
                half_patch_list[i][j] = std::make_pair(
                    patch_of_face[half_faces_list[i][j].first], half_faces_list[i][j].second);
            }
        }
        timings.push_back(timer.toc());
    }

    // group patches into shells and components
    // each shell is represented as a list of half-patch indices
    // each component is represented as a list of patch indices
    std::vector<std::vector<size_t>> shells;
    std::vector<size_t> shell_of_half_patch;
    std::vector<std::vector<size_t>> components;
    std::vector<size_t> component_of_patch;
    {
        timing_labels.emplace_back("shells and components");
        ScopedTimer<> timer("shells and components");
        compute_shells_and_components(patches.size(),
            half_patch_list,
            shells,
            shell_of_half_patch,
            components,
            component_of_patch);
        timings.push_back(timer.toc());
    }
        std::cout << "num shells = " << shells.size() << std::endl;
        std::cout << "num components = " << components.size() << std::endl;
//
    // resolve nesting order, compute material cells
    // a material cell is represented by a list of bounding shells
    std::vector<std::vector<size_t>> material_cells;
    std::vector<size_t> next_vert;
    std::vector<size_t> extremal_edge_of_component;
    {
        //        timing_labels.emplace_back("arrangement cells");
        ScopedTimer<> timer("material cells");
        if (components.size() < 2) { // no nesting problem, each shell is an arrangement cell
            material_cells.reserve(shells.size());
            for (size_t i = 0; i < shells.size(); ++i) {
                material_cells.emplace_back(1);
                material_cells.back()[0] = i;
            }
        } else { // resolve nesting order
            if (use_topo_ray_shooting) {
                // map: tet vert index --> index of next vert (with smaller (x,y,z))
                //            std::vector<size_t> next_vert;
                {
                    timing_labels.emplace_back("matCells(build next_vert)");
                    ScopedTimer<> timer("material cells: find next vert for each tet vertex");
                    build_next_vert(pts, tets, next_vert);
                    timings.push_back(timer.toc());
                }

                // find extremal edge for each component
                // extremal edge of component i is stored at position [2*i], [2*i+1]
                //            std::vector<size_t> extremal_edge_of_component;
                // store an MI-vert index on edge (v, v_next), None means there is no such MI-vert
                std::vector<size_t> MI_vert_on_v_v_next;
                // map: (tet_id, tet_face_id) --> iso_face_id
                absl::flat_hash_map<std::pair<size_t, size_t>, size_t> MI_face_id_of_tet_face;
                // map: (tet_id, tet_vert_id) --> (iso_vert_id, component_id)
                absl::flat_hash_map<std::pair<size_t, size_t>, std::pair<size_t, size_t>>
                    MI_vId_compId_of_tet_vert;
                {
                    timing_labels.emplace_back("matCells(find extremal edges)");
                    ScopedTimer<> timer("material cells: find extremal edge for components");
                    extremal_edge_of_component.resize(
                        2 * components.size(), MaterialInterface<3>::None);
                    MI_vert_on_v_v_next.resize(n_pts, MaterialInterface<3>::None);
                    MI_face_id_of_tet_face.reserve(MI_faces.size());
                    MI_vId_compId_of_tet_vert.reserve(MI_faces.size() / 2);
                    //
                    std::vector<bool> is_MI_vert_visited(MI_verts.size(), false);
                    for (size_t i = 0; i < patches.size(); ++i) {
                        size_t component_id = component_of_patch[i];
                        auto& u1 = extremal_edge_of_component[2 * component_id];
                        auto& u2 = extremal_edge_of_component[2 * component_id + 1];
                        for (auto fId : patches[i]) {
                            for (const auto& tet_face : MI_faces[fId].tet_face_indices) {
                                MI_face_id_of_tet_face.try_emplace(tet_face, fId);
                            }
                            for (auto vId : MI_faces[fId].vert_indices) {
                                if (!is_MI_vert_visited[vId]) {
                                    is_MI_vert_visited[vId] = true;
                                    const auto& vert = MI_verts[vId];
                                    if (vert.simplex_size == 2) { // edge MI-vertex
                                        auto v1 = vert.simplex_vert_indices[0];
                                        auto v2 = vert.simplex_vert_indices[1];
                                        if (next_vert[v1] == v2) { // on tree edge v1 -> v2
                                            // update extremal edge
                                            if (u1 == MaterialInterface<3>::None) {
                                                u1 = v1;
                                                u2 = v2;
                                            } else {
                                                if (v2 == u2) {
                                                    if (point_xyz_less(pts[v1], pts[u1])) {
                                                        u1 = v1;
                                                    }
                                                } else if (point_xyz_less(pts[v2], pts[u2])) {
                                                    u1 = v1;
                                                    u2 = v2;
                                                }
                                            }
                                            // record an iso-vert on edge v1 -> v2
                                            MI_vert_on_v_v_next[v1] = vId;
                                            // fill map
                                            MI_vId_compId_of_tet_vert.try_emplace(
                                                std::make_pair(vert.tet_index, vert.tet_vert_index),
                                                std::make_pair(vId, component_id));
                                        } else if (next_vert[v2] == v1) { // on tree edge v2 -> v1
                                            // update extremal edge
                                            if (u1 == MaterialInterface<3>::None) {
                                                u1 = v2;
                                                u2 = v1;
                                            } else {
                                                if (v1 == u2) {
                                                    if (point_xyz_less(pts[v2], pts[u1])) {
                                                        u1 = v2;
                                                    }
                                                } else if (point_xyz_less(pts[v1], pts[u2])) {
                                                    u1 = v2;
                                                    u2 = v1;
                                                }
                                            }
                                            // record an iso-vert on v2 -> v1
                                            MI_vert_on_v_v_next[v2] = vId;
                                            // fill map
                                            MI_vId_compId_of_tet_vert.try_emplace(
                                                std::make_pair(vert.tet_index, vert.tet_vert_index),
                                                std::make_pair(vId, component_id));
                                        }
                                    }
                                }
                            }
                        }
                    }
                    timings.push_back(timer.toc());
                }

                // topological ray shooting
                std::vector<std::pair<size_t, size_t>> shell_links;
                {
                    timing_labels.emplace_back("matCells(ray shooting)");
                    ScopedTimer<> timer("material cells: topo ray shooting");
                    shell_links.reserve(components.size());
                    std::vector<size_t> sorted_vert_indices_on_edge;
                    sorted_vert_indices_on_edge.reserve(3);
                    for (size_t i = 0; i < components.size(); ++i) {
                        // extremal edge: v1 -> v2
                        auto extreme_v1 = extremal_edge_of_component[2 * i];
                        auto extreme_v2 = extremal_edge_of_component[2 * i + 1];
                        auto MI_vId = MI_vert_on_v_v_next[extreme_v1];
                        auto tetId = MI_verts[MI_vId].tet_index;
                        const auto& tet_cut_result = cut_results[cut_result_index[tetId]];
                        // get local index of v1 and v2 in the tet
                        size_t local_v1, local_v2;
                        for (size_t j = 0; j < 4; ++j) {
                            if (tets[tetId][j] == extreme_v1) {
                                local_v1 = j;
                            } else if (tets[tetId][j] == extreme_v2) {
                                local_v2 = j;
                            }
                        }
                        //                    std::cout << "local_v1 = " << local_v1 << std::endl;
                        //                    std::cout << "local_v2 = " << local_v2 << std::endl;
                        // get an ordered list of vertices on edge v1 -> v2
                        sorted_vert_indices_on_edge.clear();
                        compute_edge_intersection_order_MI(
                            tet_cut_result, local_v1, local_v2, sorted_vert_indices_on_edge);
                        //
                        //                    std::cout << "sorted_vert_indices_on_edge" << std::endl; for (auto vId : sorted_vert_indices_on_edge) {
                        //                        std::cout << vId << ", ";
                        //                    }
                        //                    std::cout << std::endl;
                        //
                        //                    std::cout << "sorted iso-verts on edge: list of (iso-vId, compId)" << std::endl; for (int j = 1; j+1 < sorted_vert_indices_on_edge.size(); ++j) {
                        //                        const auto& iso_vId_compId =
                        //                            iso_vId_compId_of_tet_vert[std::make_pair(tetId, sorted_vert_indices_on_edge[j])];
                        //                        std::cout << "(" << iso_vId_compId.first << "," << iso_vId_compId.second << ") ";
                        //                    }
                        //                    std::cout << std::endl;
                        //
                        // find the vertex v_start on v1->v2
                        // 1. on current component
                        // 2. nearest to v2
                        size_t j_start;
                        for (size_t j = 0; j + 1 < sorted_vert_indices_on_edge.size(); ++j) {
                            const auto& MI_vId_compId = MI_vId_compId_of_tet_vert[std::make_pair(
                                tetId, sorted_vert_indices_on_edge[j])];
                            if (MI_vId_compId.second == i) {
                                j_start = j;
                            }
                        }
                        if (j_start + 2 < sorted_vert_indices_on_edge.size()) {
                            // there is a vert from another component between v_start -> v2
                            std::pair<size_t, int> face_orient1, face_orient2;
                            compute_passing_face_pair_MI(tet_cut_result,
                                sorted_vert_indices_on_edge[j_start],
                                sorted_vert_indices_on_edge[j_start + 1],
                                face_orient1,
                                face_orient2);
                            size_t MI_fId1 =
                                MI_face_id_of_tet_face[std::make_pair(tetId, face_orient1.first)];
                            size_t MI_fId2 =
                                MI_face_id_of_tet_face[std::make_pair(tetId, face_orient2.first)];
                            size_t shell1 =
                                (face_orient1.second == 1)
                                    ? shell_of_half_patch[2 * patch_of_face[MI_fId1]]
                                    : shell_of_half_patch[2 * patch_of_face[MI_fId1] + 1];
                            size_t shell2 =
                                (face_orient2.second == 1)
                                    ? shell_of_half_patch[2 * patch_of_face[MI_fId2]]
                                    : shell_of_half_patch[2 * patch_of_face[MI_fId2] + 1];
                            // link shell1 with shell2
                            shell_links.emplace_back(shell1, shell2);
                        } else {
                            // there is no vert between v_start -> v2
                            std::pair<size_t, int> face_orient;
                            compute_passing_face_MI(tet_cut_result,
                                sorted_vert_indices_on_edge[j_start],
                                sorted_vert_indices_on_edge.back(),
                                face_orient);
                            size_t MI_fId =
                                MI_face_id_of_tet_face[std::make_pair(tetId, face_orient.first)];
                            size_t shell_start =
                                (face_orient.second == 1)
                                    ? shell_of_half_patch[2 * patch_of_face[MI_fId]]
                                    : shell_of_half_patch[2 * patch_of_face[MI_fId] + 1];
                            // follow the ray till another iso-vertex or the sink
                            auto v_curr = extreme_v2;
                            while (next_vert[v_curr] != MaterialInterface<3>::None &&
                                   MI_vert_on_v_v_next[v_curr] == MaterialInterface<3>::None) {
                                v_curr = next_vert[v_curr];
                            }
                            if (MI_vert_on_v_v_next[v_curr] != MaterialInterface<3>::None) {
                                // reached iso-vert at end of the ray
                                auto MI_vId_end = MI_vert_on_v_v_next[v_curr];
                                auto end_tetId = MI_verts[MI_vId_end].tet_index;
                                const auto& end_tet_cut_result =
                                    cut_results[cut_result_index[end_tetId]];
                                auto v_next = next_vert[v_curr];
                                // find local vertex indices in the end tetrahedron
                                for (size_t j = 0; j < 4; ++j) {
                                    if (tets[end_tetId][j] == v_curr) {
                                        local_v1 = j;
                                    } else if (tets[end_tetId][j] == v_next) {
                                        local_v2 = j;
                                    }
                                }
                                // get an ordered list of vertices on edge v_curr -> v_next
                                sorted_vert_indices_on_edge.clear();
                                compute_edge_intersection_order_MI(end_tet_cut_result,
                                    local_v1,
                                    local_v2,
                                    sorted_vert_indices_on_edge);
                                // find the end shell
                                compute_passing_face_MI(end_tet_cut_result,
                                    sorted_vert_indices_on_edge[1],
                                    sorted_vert_indices_on_edge.front(),
                                    face_orient);
                                MI_fId = MI_face_id_of_tet_face[std::make_pair(
                                    end_tetId, face_orient.first)];
                                size_t shell_end =
                                    (face_orient.second == 1)
                                        ? shell_of_half_patch[2 * patch_of_face[MI_fId]]
                                        : shell_of_half_patch[2 * patch_of_face[MI_fId] + 1];
                                // link start shell with end shell
                                shell_links.emplace_back(shell_start, shell_end);
                            } else {
                                // next_vert[v_curr] is None, v_curr is the sink vertex
                                // link shell_start with the sink
                                shell_links.emplace_back(shell_start, MaterialInterface<3>::None);
                            }
                        }
                    }
                    timings.push_back(timer.toc());
                }
                //            std::cout << "shell_links = " << std::endl;
                //            for(const auto& link : shell_links) {
                //                std::cout << "(" << link.first << "," << link.second << ") ";
                //            }
                //            std::cout << std::endl;

                // group shells into arrangement cells
                {
                    timing_labels.emplace_back("matCells(group shells into matCells)");
                    ScopedTimer<> timer("material cells: group shells into material cells");
                    compute_arrangement_cells(shells.size(), shell_links, material_cells);
                    timings.push_back(timer.toc());
                }
                //            std::cout << "arrangement cells = " << std::endl;
                //            for (const auto& arr_cell : material_cells) {
                //                std::cout << "(";
                //                for (auto s : arr_cell) {
                //                    std::cout << s << ", ";
                //                }
                //                std::cout << ")" << std::endl;
                //            }
            } else { // group simplicial cells into material cells
                // ------------------- build adjacency graph of simplicial cells ---------------
                // pair (tetId, tet_cell_Id)
                std::vector<std::pair<size_t, size_t>> tet_cell_of_simp_cell;
                // simplicial half-face info
                // if the face is on material interface, let s be its shell index, record (-s-1)
                // if the face is not on material interface, record the simp_cell index on the opposite side
                std::vector<long long> simp_half_face_info;
                std::vector<size_t> simp_hFace_start_index;
                {
                    timing_labels.emplace_back("matCells(build simpCell graph)");
                    ScopedTimer<> timer("material cells: build simplicial cells graph");
                    size_t num_simplicial_cells = 0;
                    size_t num_simp_half_faces = 0;
                    for (size_t i = 0; i < tets.size(); i++) {
                        if (cut_result_index[i] != MaterialInterface<3>::None) {
                            const auto& cells = cut_results[cut_result_index[i]].cells;
                            num_simplicial_cells += cells.size();
                            for (const auto& cell : cells) {
                                num_simp_half_faces += cell.faces.size();
                            }
                        } else { // empty tet
                            ++num_simplicial_cells;
                            num_simp_half_faces += 4;
                        }
                    }
                    tet_cell_of_simp_cell.reserve(num_simplicial_cells);
                    simp_half_face_info.reserve(num_simp_half_faces);
                    simp_hFace_start_index.reserve(num_simplicial_cells + 1);
                    simp_hFace_start_index.push_back(0);
                    // hash table: face key --> (simplicial cell index, face index in simplicial
                    // cell)
                    absl::flat_hash_map<std::array<long long, 3>, std::pair<size_t, size_t>>
                        incident_cell_of_face;
                    //                    incident_cell_of_face.reserve();
                    std::vector<long long> face_verts;
                    face_verts.reserve(4);
                    std::array<long long, 3> key;
                    std::vector<std::array<long long, 3>> four_face_verts(4);
                    for (size_t i = 0; i < tets.size(); i++) {
                        if (cut_result_index[i] != MaterialInterface<3>::None) {
                            size_t MI_fId_start_index = MI_fId_start_index_of_tet[i];
                            size_t global_vId_start_index = global_vId_start_index_of_tet[i];
                            const auto& matInterface = cut_results[cut_result_index[i]];
                            const auto& faces = matInterface.faces;
                            const auto& cells = matInterface.cells;
                            for (size_t j = 0; j < cells.size(); j++) {
                                // create a new simplicial cell
                                const auto& cell = cells[j];
                                size_t cur_simp_cell_id = tet_cell_of_simp_cell.size();
                                tet_cell_of_simp_cell.emplace_back(i, j);
                                for (size_t k = 0; k < cell.faces.size(); k++) {
                                    size_t fId = cell.faces[k];
                                    size_t MI_fId = MI_fId_of_tet_face[MI_fId_start_index + fId];
                                    if (MI_fId !=
                                        MaterialInterface<3>::None) { // face k is on iso-surface
                                        size_t patch_id = patch_of_face[MI_fId];
                                        size_t half_patch_id = (faces[fId].positive_material_label == cell.material_label)
                                                                   ? (2 * patch_id)
                                                                   : (2 * patch_id + 1);
                                        size_t shell_id = shell_of_half_patch[half_patch_id];
                                        simp_half_face_info.push_back(-shell_id - 1);
                                    } else { // face k is not on material interface
                                        face_verts.clear();
                                        const auto& faceVertices = faces[fId].vertices;
                                        // convert from local vert index to global vert index
                                        for (unsigned long vId : faceVertices) {
                                            face_verts.push_back(
                                                global_vId_of_tet_vert[global_vId_start_index + vId]);
                                        }
                                        //
                                        compute_iso_face_key(face_verts, key);
                                        auto iter_inserted = incident_cell_of_face.try_emplace(
                                            key, std::make_pair(cur_simp_cell_id, k));
                                        if (!iter_inserted
                                                 .second) { // same face has been inserted before
                                            size_t opposite_simp_cell_id =
                                                iter_inserted.first->second.first;
                                            size_t opposite_cell_face_id =
                                                iter_inserted.first->second.second;
                                            // make neighbor
                                            simp_half_face_info.push_back(opposite_simp_cell_id);
                                            simp_half_face_info
                                                [simp_hFace_start_index[opposite_simp_cell_id] +
                                                    opposite_cell_face_id] = cur_simp_cell_id;
                                            // delete face in hash table since a face can only be
                                            // shared by two cells
                                            incident_cell_of_face.erase(iter_inserted.first);
                                        } else { // the face is inserted for the first time
                                            simp_half_face_info.push_back(
                                                std::numeric_limits<long long>::max());
                                        }
                                    }
                                }
                                simp_hFace_start_index.push_back(simp_half_face_info.size());
                            }
                        } else { // tet i has no material interface
                            // create a new simplicial cell
                            size_t cur_simp_cell_id = tet_cell_of_simp_cell.size();
                            tet_cell_of_simp_cell.emplace_back(i, 0);
                            // global index of four tet vertices
                            long long global_vId0 = -tets[i][0] - 1;
                            long long global_vId1 = -tets[i][1] - 1;
                            long long global_vId2 = -tets[i][2] - 1;
                            long long global_vId3 = -tets[i][3] - 1;
                            // face 0
                            four_face_verts[0] = {global_vId1, global_vId2, global_vId3};
                            std::sort(four_face_verts[0].begin(), four_face_verts[0].end());
                            // face 1
                            four_face_verts[1] = {global_vId2, global_vId3, global_vId0};
                            std::sort(four_face_verts[1].begin(), four_face_verts[1].end());
                            // face 2
                            four_face_verts[2] = {global_vId3, global_vId0, global_vId1};
                            std::sort(four_face_verts[2].begin(), four_face_verts[2].end());
                            // face 3
                            four_face_verts[3] = {global_vId0, global_vId1, global_vId2};
                            std::sort(four_face_verts[3].begin(), four_face_verts[3].end());
                            //
                            for (size_t j = 0; j < 4; j++) {
                                auto iter_inserted = incident_cell_of_face.try_emplace(
                                    four_face_verts[j], std::make_pair(cur_simp_cell_id, j));
                                if (!iter_inserted.second) { // same face has been inserted before
                                    size_t opposite_simp_cell_id =
                                        iter_inserted.first->second.first;
                                    size_t opposite_cell_face_id =
                                        iter_inserted.first->second.second;
                                    // make neighbor
                                    simp_half_face_info.push_back(opposite_simp_cell_id);
                                    simp_half_face_info
                                        [simp_hFace_start_index[opposite_simp_cell_id] +
                                            opposite_cell_face_id] = cur_simp_cell_id;
                                    // delete face in hash table since a face can only be shared by
                                    // two cells
                                    incident_cell_of_face.erase(iter_inserted.first);
                                } else { // face inserted for the first time
                                    simp_half_face_info.push_back(
                                        std::numeric_limits<long long>::max());
                                }
                            }
                            simp_hFace_start_index.push_back(simp_half_face_info.size());
                        }
                    }
                    timings.push_back(timer.toc());
                }
                {
                    // ------------------- group simplicial cells into arrangement cells
                    // ---------------
                    timing_labels.emplace_back("matCells(group simpCells into matCells)");
                    ScopedTimer<> timer("material cells: group simplicial cells");
                    size_t num_simp_cells = tet_cell_of_simp_cell.size();
                    std::vector<bool> visited_simp_cell(num_simp_cells, false);
                    std::vector<absl::flat_hash_set<size_t>> material_cell_incident_shells;
                    for (size_t i = 0; i < num_simp_cells; i++) {
                        if (!visited_simp_cell[i]) {
                            // new material cell
                            material_cell_incident_shells.emplace_back();
                            auto& incident_shells = material_cell_incident_shells.back();
                            //
                            std::queue<size_t> Q;
                            Q.push(i);
                            visited_simp_cell[i] = true;
                            while (!Q.empty()) {
                                size_t simp_cell_id = Q.front();
                                Q.pop();
                                size_t start_index = simp_hFace_start_index[simp_cell_id];
                                size_t face_info_size =
                                    simp_hFace_start_index[simp_cell_id + 1] - start_index;
                                for (size_t j = 0; j < face_info_size; j++) {
                                    if (simp_half_face_info[start_index + j] < 0) {
                                        // the half-face is on material interface
                                        size_t shell_id = -simp_half_face_info[start_index + j] - 1;
                                        incident_shells.insert(shell_id);
                                    } else {
                                        // the half-face is not on material interface
                                        long long opposite_simp_cell_id =
                                            simp_half_face_info[start_index + j];
                                        if ((opposite_simp_cell_id !=
                                                std::numeric_limits<long long>::max()) &&
                                            !visited_simp_cell[opposite_simp_cell_id]) {
                                            // the opposite simplicial cell exists and is not
                                            // visited
                                            Q.push(opposite_simp_cell_id);
                                            visited_simp_cell[opposite_simp_cell_id] = true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // convert vector-of-set to vector-of-vector
                    material_cells.reserve(material_cell_incident_shells.size());
                    for (auto & material_cell_incident_shell : material_cell_incident_shells) {
                        material_cells.emplace_back(material_cell_incident_shell.begin(),
                            material_cell_incident_shell.end());
                    }
                    timings.push_back(timer.toc());
                }
            }
        }
        timings.push_back(timer.toc());
        timing_labels.emplace_back("material cells");
    }
    std::cout << "num_cells = " << material_cells.size() << std::endl;
    if (components.size() > 1) {
        timing_labels.emplace_back("matCells(other)");
        size_t num_timings = timings.size();
        if (use_topo_ray_shooting) {
            timings.push_back(timings[num_timings - 1] - timings[num_timings - 2] -
                              timings[num_timings - 3] - timings[num_timings - 4] -
                              timings[num_timings - 5] - timings[num_timings - 6]);
        } else {
            // baseline: group simplicial cells into material cells
            timings.push_back(
                timings[num_timings - 1] - timings[num_timings - 2] - timings[num_timings - 3]);
        }
    }


    // test: export mesh, patches, chains
    if (!args.timing_only) {
        save_result_MI(output_dir + "/MI_mesh.json",
            MI_pts,
            MI_faces,
            patches,
            MI_edges,
            chains,
            non_manifold_edges_of_vert,
            half_patch_list,
            shells,
            components,
            material_cells);
        save_result_msh(output_dir + "/MI_mesh",
            MI_pts,
            MI_faces,
            patches,
            MI_edges,
            chains,
            non_manifold_edges_of_vert,
            half_patch_list,
            shells,
            components,
            material_cells);
        //
//        if (components.size() > 1) {
//            save_nesting_data(output_dir + "/nesting_data.json",
//                next_vert,
//                extremal_edge_of_component);
//        }
    }
    // test: export timings
    save_timings(output_dir + "/timings.json", timing_labels, timings);

    return 0;
}
