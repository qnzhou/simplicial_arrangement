#include <simplicial_arrangement/lookup_table.h>
#include <simplicial_arrangement/simplicial_arrangement.h>

#include <chrono>
#include <iostream>
#include <exception>
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
    CLI::App app{"Implicit Arrangement Command Line"};
    app.add_option("config_file", args.config_file, "Configuration file")->required();
    app.add_option("-T,--timing-only", args.timing_only, "Record timing without output result");
    app.add_option("-R,--robust-test",args.robust_test, "Perform robustness test");
    CLI11_PARSE(app, argc, argv);

    // parse configure file
    std::string tet_mesh_file;
//    std::string sphere_file;
    std::string func_file;
    std::string output_dir;
    bool use_lookup = true;
    bool use_2func_lookup = true;
    bool use_topo_ray_shooting = true;
    bool use_bbox = true;
    std::array<double, 3> bbox_min, bbox_max;
    parse_config_file(args.config_file,
        tet_mesh_file,
        func_file,
        output_dir,
        use_lookup,
        use_2func_lookup,
        use_topo_ray_shooting,
        use_bbox,
        bbox_min,
        bbox_max);
    //    std::string config_path = args.config_file.substr(0, args.config_file.find_last_of('/'));
    //    std::cout << "config path: " << config_path << std::endl;
    //    tet_mesh_file = config_path + "/" + tet_mesh_file;
    //    sphere_file = config_path + "/" + sphere_file;
    if (use_lookup) {
        // load lookup table
        std::cout << "load table ..." << std::endl;
        bool loaded = load_lookup_table();
        if (loaded) {
            std::cout << "loading finished." << std::endl;
        } else {
            std::cout << "loading failed." << std::endl;
            return -1;
        }
    } else {
        use_2func_lookup = false;
    }

    // record timings
    std::vector<std::string> timing_labels;
    std::vector<double> timings;

    // record stats
    std::vector<std::string> stats_labels;
    std::vector<size_t> stats;

    // load tet mesh
    std::vector<std::array<double, 3>> pts;
    std::vector<std::array<size_t, 4>> tets;
    load_tet_mesh(tet_mesh_file, pts, tets);
    size_t n_tets = tets.size();
    size_t n_pts = pts.size();
    std::cout << "tet mesh: " << n_pts << " verts, " << n_tets << " tets." << std::endl;
    stats_labels.emplace_back("num_pts");
    stats.push_back(n_pts);
    stats_labels.emplace_back("num_tets");
    stats.push_back(n_tets);


    // load implicit function values, or evaluate
//    std::vector<Sphere> spheres;
//    load_spheres(sphere_file, spheres);
//    size_t n_func = spheres.size();
//
//
//
//    {
//        timing_labels.emplace_back("func values");
//        ScopedTimer<> timer("func values");
//        funcVals.resize(n_pts, n_func);
//        for (Eigen::Index i = 0; i < n_pts; ++i) {
//            const auto& p = pts[i];
//            for (Eigen::Index j = 0; j < n_func; ++j) {
//                funcVals(i,j) = sphere_function(spheres[j].first, spheres[j].second, p);
//            }
//        }
//        timings.push_back(timer.toc());
//    }

//    if (args.robust_test) {
    if (false) {
        // load robustness test spheres
        std::vector<Sphere> spheres;
        load_spheres(func_file, spheres);

        size_t n_func = 4;
        size_t n_test = spheres.size() / n_func;
//        size_t n_test = 1;
        //
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> funcVals(n_pts, n_func);
        Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> funcSigns(n_pts, n_func);
        std::vector<bool> is_degenerate_vertex;
        is_degenerate_vertex.reserve(n_pts);
        //
        std::vector<size_t> func_in_tet;
        func_in_tet.reserve(n_tets);
        std::vector<size_t> start_index_of_tet;
        start_index_of_tet.reserve(n_tets+1);
        //
        Arrangement<3> cut_result, cut_result2;
        std::vector<Plane<double, 3>> planes;
        std::vector<Plane<double, 3>> planes_reverse; // planes in reverse order
        planes.reserve(3);
        planes_reverse.reserve(3);
        //
        size_t type1_count = 0;  // normal and reverse order run, but are different
        size_t type2_count = 0;  // normal order crash
        size_t type3_count = 0;  // normal order run, but reverse order crash
        //
        for (size_t iter = 0; iter < n_test; ++iter) {
            std::cout << "-------- test " << iter << " -----------" << std::endl;

            // load implicit functions
            for (Eigen::Index i = 0; i < n_pts; ++i) {
                const auto& p = pts[i];
                for (Eigen::Index j = 0; j < n_func; ++j) {
                    funcVals(i, j) = compute_sphere_distance(spheres[iter* n_func + j].first,
                        spheres[iter* n_func + j].second, p);
                }
            }
//            load_functions(func_file, pts, funcVals);
            std::cout << "funcVals(0,0) = " << funcVals(0,0) << std::endl;
            // function signs at vertices
            bool found_degenerate_vertex = false;
            size_t num_degenerate_vertex = 0;
            {
                timing_labels.emplace_back("func signs");
                ScopedTimer<> timer("func signs");
                is_degenerate_vertex.clear();
                is_degenerate_vertex.resize(n_pts, false);
                for (Eigen::Index i = 0; i < n_pts; i++) {
                    for (Eigen::Index j = 0; j < n_func; j++) {
                        funcSigns(i, j) = sign(funcVals(i, j));
                        if (funcSigns(i, j) == 0) {
                            is_degenerate_vertex[i] = true;
                            found_degenerate_vertex = true;
                            num_degenerate_vertex++;
                        }
                    }
                }
                timings.push_back(timer.toc());
            }
            std::cout << "num_degenerate_vertex = " << num_degenerate_vertex << std::endl;
            // filter
            size_t num_intersecting_tet = 0;
            {
                timing_labels.emplace_back("filter");
                ScopedTimer<> timer("filter(CRS vector)");
//                func_in_tet.reserve(n_tets);
//                start_index_of_tet.reserve(n_tets + 1);
                func_in_tet.clear();
                start_index_of_tet.clear();
                start_index_of_tet.push_back(0);
                int pos_count;
                int neg_count;
                for (Eigen::Index i = 0; i < n_tets; i++) {
                    for (Eigen::Index j = 0; j < n_func; j++) {
                        pos_count = 0;
                        neg_count = 0;
                        for (size_t& vId : tets[i]) {
                            if (funcSigns( vId,j) == 1) {
                                pos_count += 1;
                            } else if (funcSigns( vId,j) == -1) {
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
            Time_duration time_1_func = Time_duration::zero();
            Time_duration time_2_func = Time_duration::zero();
            Time_duration time_more_func = Time_duration::zero();
            size_t num_1_func = 0;
            size_t num_2_func = 0;
            size_t num_more_func = 0;
            //
            bool is_type2 = false;
            bool is_type1 = false;
            bool is_type3 = false;
            //
            {
                timing_labels.emplace_back("simp_arr(other)");
                ScopedTimer<> timer("simp_arr");
                size_t start_index;
                size_t num_func;
                for (size_t i = 0; i < tets.size(); i++) {
                    start_index = start_index_of_tet[i];
                    num_func = start_index_of_tet[i + 1] - start_index;
                    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
                    //
                    size_t v1 = tets[i][0];
                    size_t v2 = tets[i][1];
                    size_t v3 = tets[i][2];
                    size_t v4 = tets[i][3];
                    planes.clear();
                    for (size_t j = 0; j < num_func; j++) {
                        size_t f_id = func_in_tet[start_index + j];
                        planes.emplace_back();
                        auto& plane = planes.back();
                        plane[0] = funcVals(v1, f_id);
                        plane[1] = funcVals(v2, f_id);
                        plane[2] = funcVals(v3, f_id);
                        plane[3] = funcVals(v4, f_id);
                    }
                    // reverse plane order
                    planes_reverse.clear();
                    for (size_t j = 0; j < num_func; ++j) {
                        planes_reverse.emplace_back(planes[num_func - 1 - j]);
                    }
                    //
                    bool crashed = false;
                    if (!use_2func_lookup && num_func == 2) {
                        disable_lookup_table();
                        try {
                            cut_result = compute_arrangement(planes);
                        } catch (std::runtime_error& e) {
                            crashed = true;
                            is_type2 = true;
                            break;
                        }
                        if (!crashed) {
                            try {
                                cut_result2 = compute_arrangement(planes_reverse);
                            } catch (std::runtime_error& e) {
                                crashed = true;
//                                is_type2 = true;
                                is_type3 = true;
//                                break;
                            }
                            if (!crashed) {
                                if (cut_result.vertices.size() !=
                                        cut_result2.vertices.size() ||
                                    cut_result.faces.size() != cut_result2.faces.size() ||
                                    cut_result.cells.size() != cut_result2.cells.size()) {
                                    // inconsistent results
                                    is_type1 = true;
//                                    break;
                                }
                            }
                        }
                        enable_lookup_table();
                    } else {
                        try {
                            cut_result =compute_arrangement(planes);
                        } catch (std::runtime_error& e) {
                            crashed = true;
                            is_type2 = true;
                            break;
                        }
                        if (!crashed) {
                            try {
                                cut_result2 = compute_arrangement(planes_reverse);
                            } catch (std::runtime_error& e) {
                                crashed = true;
//                                is_type2 = true;
                                is_type3 = true;
//                                break;
                            }
                            if (!crashed) {
                                if (cut_result.vertices.size() !=
                                        cut_result2.vertices.size() ||
                                    cut_result.faces.size() != cut_result2.faces.size() ||
                                    cut_result.cells.size() != cut_result2.cells.size()) {
                                    // inconsistent results
                                    is_type1 = true;
                                }
                            }
                        }
                    }
                    //
                    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
                    switch (num_func) {
                    case 1:
                        time_1_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                        ++num_1_func;
                        break;
                    case 2:
                        time_2_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                        ++num_2_func;
                        break;
                    default:
                        time_more_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                        ++num_more_func;
                        break;
                    }
                }
                timings.push_back(timer.toc() - time_1_func.count() - time_2_func.count() -
                                  time_more_func.count());
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
    if (load_functions(func_file, pts, funcVals)) {
        std::cout << "function loading finished." << std::endl;
    } else {
        std::cout << "function loading failed." << std::endl;
        return -2;
    }
    size_t n_func = funcVals.cols();


    // function signs at vertices
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> funcSigns;
    std::vector<bool> is_degenerate_vertex;
    bool found_degenerate_vertex = false;
    size_t num_degenerate_vertex = 0;
    {
        timing_labels.emplace_back("func signs");
        ScopedTimer<> timer("func signs");
        is_degenerate_vertex.resize(n_pts, false);
        funcSigns.resize(n_pts, n_func);
        for (Eigen::Index i = 0; i < n_pts; i++) {
            for (Eigen::Index j = 0; j < n_func; j++) {
                funcSigns(i, j) = sign(funcVals(i, j));
                if (funcSigns(i, j) == 0) {
                    is_degenerate_vertex[i] = true;
                    found_degenerate_vertex = true;
                    num_degenerate_vertex++;
                }
            }
        }
        timings.push_back(timer.toc());
    }
    std::cout << "num_degenerate_vertex = " << num_degenerate_vertex << std::endl;
    stats_labels.emplace_back("num_degenerate_vertex");
    stats.push_back(num_degenerate_vertex);

    // filter
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
                    if (funcSigns( vId,j) == 1) {
                        pos_count += 1;
                    } else if (funcSigns( vId,j) == -1) {
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
    stats_labels.emplace_back("num_intersecting_tet");
    stats.push_back(num_intersecting_tet);


    // compute arrangement in each tet (iterative plane cut)
    std::vector<Arrangement<3>> cut_results;
    std::vector<size_t> cut_result_index;
    //
    Time_duration time_1_func = Time_duration::zero();
    Time_duration time_2_func = Time_duration::zero();
    Time_duration time_more_func = Time_duration::zero();
    size_t num_1_func = 0;
    size_t num_2_func = 0;
    size_t num_more_func = 0;
    //
    bool is_type2 = false;
    bool is_type3 = false;
    bool is_type1 = false;
    //
    {
        timing_labels.emplace_back("simp_arr(other)");
        ScopedTimer<> timer("simp_arr");
        if (args.robust_test) {
            Arrangement<3> cut_result, cut_result2;
            size_t start_index;
            size_t num_func;
            std::vector<Plane<double, 3>> planes;
            std::vector<Plane<double, 3>> planes_reverse; // planes in reverse order
            planes.reserve(3);
            planes_reverse.reserve(3);
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
                planes.clear();
                for (size_t j = 0; j < num_func; j++) {
                    size_t f_id = func_in_tet[start_index + j];
                    planes.emplace_back();
                    auto& plane = planes.back();
                    plane[0] = funcVals(v1, f_id);
                    plane[1] = funcVals(v2, f_id);
                    plane[2] = funcVals(v3, f_id);
                    plane[3] = funcVals(v4, f_id);
                }
                // reverse plane order
                planes_reverse.clear();
                for (size_t j = 0; j < num_func; ++j) {
                    planes_reverse.emplace_back(planes[num_func -1 -j]);
                }
                //
                bool crashed = false;
                if (!use_2func_lookup && num_func == 2) {
                    cut_result_index.push_back(cut_results.size());
                    disable_lookup_table();
                    try {
                        cut_result = compute_arrangement(planes);
                    } catch (std::runtime_error& e) {
                        crashed = true;
                        is_type2 = true;
                        break;
                    }
                    if (!crashed) {
                        try {
                            cut_result2 = compute_arrangement(planes_reverse);
                        } catch (std::runtime_error& e) {
                            crashed = true;
                            is_type3 = true;
                        }
                        if (!crashed) {
                            if (cut_result.vertices.size() != cut_result2.vertices.size() ||
                                cut_result.faces.size() != cut_result2.faces.size() ||
                                cut_result.cells.size() != cut_result2.cells.size()) {
                                // inconsistent results
                                is_type1 = true;
                            }
                        }
                    }
                    enable_lookup_table();
                } else {
                    try {
                        cut_result = compute_arrangement(planes);
                    } catch (std::runtime_error& e) {
                        crashed = true;
                        is_type2 = true;
                        break;
                    }
                    if (!crashed) {
                        try {
                            cut_result2 = compute_arrangement(planes_reverse);
                        } catch (std::runtime_error& e) {
                            crashed = true;
                            is_type3 = true;
                        }
                        if (!crashed) {
                            if (cut_result.vertices.size() != cut_result2.vertices.size() ||
                                cut_result.faces.size() != cut_result2.faces.size() ||
                                cut_result.cells.size() != cut_result2.cells.size()) {
                                // inconsistent results
                                is_type1 = true;
                            }
                        }
                    }
                }
                //
                std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
                switch (num_func) {
                case 1:
                    time_1_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                    ++num_1_func;
                    break;
                case 2:
                    time_2_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                    ++num_2_func;
                    break;
                default:
                    time_more_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                    ++num_more_func;
                    break;
                }
            }
        } else { // not performing robustness test
            cut_results.reserve(num_intersecting_tet);
            cut_result_index.reserve(n_tets);
            size_t start_index;
            size_t num_func;
            std::vector<Plane<double, 3>> planes;
            planes.reserve(3);
            try {
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
                    planes.clear();
                    for (size_t j = 0; j < num_func; j++) {
                        size_t f_id = func_in_tet[start_index + j];
                        planes.emplace_back();
                        auto& plane = planes.back();
                        plane[0] = funcVals(v1, f_id);
                        plane[1] = funcVals(v2, f_id);
                        plane[2] = funcVals(v3, f_id);
                        plane[3] = funcVals(v4, f_id);
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
                    case 1:
                        time_1_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                        ++num_1_func;
                        break;
                    case 2:
                        time_2_func += std::chrono::duration_cast<Time_duration>(t2 - t1);
                        ++num_2_func;
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
    //
    stats_labels.emplace_back("num_1_func");
    stats.push_back(num_1_func);
    stats_labels.emplace_back("num_2_func");
    stats.push_back(num_2_func);
    stats_labels.emplace_back("num_more_func");
    stats.push_back(num_more_func);

    if (args.robust_test) {
        if (is_type2) {
            std::cout << "type 2 failure." << std::endl;
        } else if (is_type3) {
            std::cout << "type 3 failure." << std::endl;
        } else if (is_type1) {
            std::cout << "type 1 failure." << std::endl;
        } else {
            std::cout << "success." << std::endl;
        }
        return 0;
    }


    // extract arrangement mesh
    std::vector<IsoVert> iso_verts;
    std::vector<PolygonFace> iso_faces;
    // the following data are only needed when we use
    // the baseline nesting algorithm
    // (group simplicial cells into arrangement cells)
    std::vector<long long> global_vId_of_tet_vert;
    std::vector<size_t> global_vId_start_index_of_tet;
    std::vector<size_t> iso_fId_of_tet_face;
    std::vector<size_t> iso_fId_start_index_of_tet;
    {
        timing_labels.emplace_back("extract mesh");
        ScopedTimer<> timer("extract mesh");
        if (use_topo_ray_shooting) {
            extract_iso_mesh_pure(num_1_func,
                num_2_func,
                num_more_func,
                cut_results,
                cut_result_index,
                func_in_tet,
                start_index_of_tet,
                tets,
                iso_verts,
                iso_faces);
        } else { // nesting algorithm: group simplicial cells into arrangement cells
            extract_iso_mesh(num_1_func,
                num_2_func,
                num_more_func,
                cut_results,
                cut_result_index,
                func_in_tet,
                start_index_of_tet,
                tets,
                iso_verts,
                iso_faces,
                global_vId_of_tet_vert,
                global_vId_start_index_of_tet,
                iso_fId_of_tet_face,
                iso_fId_start_index_of_tet);
        }
        timings.push_back(timer.toc());
    }
    std::cout << "num iso-vertices = " << iso_verts.size() << std::endl;
    std::cout << "num iso-faces = " << iso_faces.size() << std::endl;
    stats_labels.emplace_back("num_iso_verts");
    stats.push_back(iso_verts.size());
    stats_labels.emplace_back("num_iso_faces");
    stats.push_back(iso_faces.size());

    // compute xyz of iso-vertices
    std::vector<std::array<double, 3>> iso_pts;
//    std::vector<std::array<long double, 3>> iso_pts_ldouble;
    {
        timing_labels.emplace_back("compute xyz");
        ScopedTimer<> timer("compute xyz");
        compute_iso_vert_xyz(iso_verts, funcVals, pts, iso_pts);
        timings.push_back(timer.toc());
    }
//    iso_pts.resize(iso_pts_ldouble.size());
//    for (size_t i = 0; i < iso_pts.size(); ++i) {
//        iso_pts[i][0] = iso_pts_ldouble[i][0];
//        iso_pts[i][1] = iso_pts_ldouble[i][1];
//        iso_pts[i][2] = iso_pts_ldouble[i][2];
//    }

    //  compute iso-edges and edge-face connectivity
    std::vector<Edge> iso_edges;
    std::vector<std::vector<size_t>> edges_of_iso_face;
    {
        timing_labels.emplace_back("isoEdge-face connectivity");
        ScopedTimer<> timer("isoEdge-face connectivity");
        compute_mesh_edges(iso_faces, edges_of_iso_face, iso_edges);
        timings.push_back(timer.toc());
    }
    std::cout << "num iso-edges = " << iso_edges.size() << std::endl;
    stats_labels.emplace_back("num_iso_edges");
    stats.push_back(iso_edges.size());


    // group iso-faces into patches
    std::vector<std::vector<size_t>> patches;
    {
        timing_labels.emplace_back("patches");
        ScopedTimer<> timer("patches");
        compute_patches(edges_of_iso_face, iso_edges, patches);
        timings.push_back(timer.toc());
    }
    std::cout << "num patches = " << patches.size() << std::endl;
    stats_labels.emplace_back("num_patches");
    stats.push_back(patches.size());

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
    std::cout << "num chains = " << chains.size() << std::endl;
    stats_labels.emplace_back("num_chains");
    stats.push_back(chains.size());


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
    stats_labels.emplace_back("num_shells");
    stats.push_back(shells.size());
    stats_labels.emplace_back("num_components");
    stats.push_back(components.size());

    // resolve nesting order, compute arrangement cells
    // an arrangement cell is represented by a list of bounding shells
    std::vector<std::vector<size_t>> arrangement_cells;
    std::vector<size_t> next_vert;
    std::vector<size_t> extremal_edge_of_component;
    {
        ScopedTimer<> timer("arrangement cells");
        if (components.size() < 2) { // no nesting problem, each shell is an arrangement cell
            arrangement_cells.reserve(shells.size());
            for (size_t i = 0; i < shells.size(); ++i) {
                arrangement_cells.emplace_back(1);
                arrangement_cells.back()[0] = i;
            }
        } else { // resolve nesting order
            if (use_topo_ray_shooting) {
                // map: tet vert index --> index of next vert (with smaller (x,y,z))
                //            std::vector<size_t> next_vert;
                {
                    timing_labels.emplace_back("arrCells(build next_vert)");
                    ScopedTimer<> timer("arrangement cells: find next vert for each tet vertex");
                    build_next_vert(pts, tets, next_vert);
                    timings.push_back(timer.toc());
                }

                //
                // find extremal edge for each component
                // extremal edge of component i is stored at position [2*i], [2*i+1]
                //            std::vector<size_t> extremal_edge_of_component;
                // store an iso-vert index on edge (v, v_next), None means there is no such iso-vert
                std::vector<size_t> iso_vert_on_v_v_next;
                // map: (tet_id, tet_face_id) --> iso_face_id
                absl::flat_hash_map<std::pair<size_t, size_t>, size_t> iso_face_id_of_tet_face;
                // map: (tet_id, tet_vert_id) --> (iso_vert_id, component_id)
                absl::flat_hash_map<std::pair<size_t, size_t>, std::pair<size_t, size_t>>
                    iso_vId_compId_of_tet_vert;
                {
                    timing_labels.emplace_back("arrCells(find extremal edges)");
                    ScopedTimer<> timer("arrangement cells: find extremal edge for components");
                    extremal_edge_of_component.resize(2 * components.size(), Arrangement<3>::None);
                    iso_vert_on_v_v_next.resize(n_pts, Arrangement<3>::None);
                    iso_face_id_of_tet_face.reserve(iso_faces.size());
                    iso_vId_compId_of_tet_vert.reserve(iso_faces.size() / 2);
                    //
                    std::vector<bool> is_iso_vert_visited(iso_verts.size(), false);
                    for (size_t i = 0; i < patches.size(); ++i) {
                        size_t component_id = component_of_patch[i];
                        auto& u1 = extremal_edge_of_component[2 * component_id];
                        auto& u2 = extremal_edge_of_component[2 * component_id + 1];
                        for (auto fId : patches[i]) {
                            for (const auto& tet_face : iso_faces[fId].tet_face_indices) {
                                iso_face_id_of_tet_face.try_emplace(tet_face, fId);
                            }
                            for (auto vId : iso_faces[fId].vert_indices) {
                                if (!is_iso_vert_visited[vId]) {
                                    is_iso_vert_visited[vId] = true;
                                    const auto& vert = iso_verts[vId];
                                    if (vert.simplex_size == 2) { // edge iso-vertex
                                        auto v1 = vert.simplex_vert_indices[0];
                                        auto v2 = vert.simplex_vert_indices[1];
                                        if (next_vert[v1] == v2) { // on tree edge v1 -> v2
                                            // update extremal edge
                                            if (u1 == Arrangement<3>::None) {
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
                                            iso_vert_on_v_v_next[v1] = vId;
                                            // fill map
                                            iso_vId_compId_of_tet_vert.try_emplace(
                                                std::make_pair(vert.tet_index, vert.tet_vert_index),
                                                std::make_pair(vId, component_id));
                                        } else if (next_vert[v2] == v1) { // on tree edge v2 -> v1
                                            // update extremal edge
                                            if (u1 == Arrangement<3>::None) {
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
                                            iso_vert_on_v_v_next[v2] = vId;
                                            // fill map
                                            iso_vId_compId_of_tet_vert.try_emplace(
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
                    timing_labels.emplace_back("arrCells(ray shooting)");
                    ScopedTimer<> timer("arrangement cells: topo ray shooting");
                    shell_links.reserve(components.size());
                    std::vector<size_t> sorted_vert_indices_on_edge;
                    sorted_vert_indices_on_edge.reserve(3);
                    for (size_t i = 0; i < components.size(); ++i) {
                        // extremal edge: v1 -> v2
                        auto extreme_v1 = extremal_edge_of_component[2 * i];
                        auto extreme_v2 = extremal_edge_of_component[2 * i + 1];
                        auto iso_vId = iso_vert_on_v_v_next[extreme_v1];
                        auto tetId = iso_verts[iso_vId].tet_index;
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
                        compute_edge_intersection_order(
                            tet_cut_result, local_v1, local_v2, sorted_vert_indices_on_edge);
                        //
                        //                    std::cout << "sorted_vert_indices_on_edge" <<
                        //                    std::endl; for (auto vId :
                        //                    sorted_vert_indices_on_edge) {
                        //                        std::cout << vId << ", ";
                        //                    }
                        //                    std::cout << std::endl;
                        //
                        //                    std::cout << "sorted iso-verts on edge: list of
                        //                    (iso-vId, compId)" << std::endl; for (int j = 1; j+1 <
                        //                    sorted_vert_indices_on_edge.size(); ++j) {
                        //                        const auto& iso_vId_compId =
                        //                            iso_vId_compId_of_tet_vert[std::make_pair(tetId,
                        //                            sorted_vert_indices_on_edge[j])];
                        //                        std::cout << "(" << iso_vId_compId.first << "," <<
                        //                        iso_vId_compId.second << ") ";
                        //                    }
                        //                    std::cout << std::endl;
                        //
                        // find the vertex v_start on v1->v2
                        // 1. on current component
                        // 2. nearest to v2
                        size_t j_start;
                        for (size_t j = 0; j + 1 < sorted_vert_indices_on_edge.size(); ++j) {
                            const auto& iso_vId_compId = iso_vId_compId_of_tet_vert[std::make_pair(
                                tetId, sorted_vert_indices_on_edge[j])];
                            if (iso_vId_compId.second == i) {
                                j_start = j;
                            }
                        }
                        if (j_start + 2 < sorted_vert_indices_on_edge.size()) {
                            // there is a vert from another component between v_start -> v2
                            std::pair<size_t, int> face_orient1, face_orient2;
                            compute_passing_face_pair(tet_cut_result,
                                sorted_vert_indices_on_edge[j_start],
                                sorted_vert_indices_on_edge[j_start + 1],
                                face_orient1,
                                face_orient2);
                            size_t iso_fId1 =
                                iso_face_id_of_tet_face[std::make_pair(tetId, face_orient1.first)];
                            size_t iso_fId2 =
                                iso_face_id_of_tet_face[std::make_pair(tetId, face_orient2.first)];
                            size_t shell1 =
                                (face_orient1.second == 1)
                                    ? shell_of_half_patch[2 * patch_of_face[iso_fId1]]
                                    : shell_of_half_patch[2 * patch_of_face[iso_fId1] + 1];
                            size_t shell2 =
                                (face_orient2.second == 1)
                                    ? shell_of_half_patch[2 * patch_of_face[iso_fId2]]
                                    : shell_of_half_patch[2 * patch_of_face[iso_fId2] + 1];
                            // link shell1 with shell2
                            shell_links.emplace_back(shell1, shell2);
                        } else {
                            // there is no vert between v_start -> v2
                            std::pair<size_t, int> face_orient;
                            compute_passing_face(tet_cut_result,
                                sorted_vert_indices_on_edge[j_start],
                                sorted_vert_indices_on_edge.back(),
                                face_orient);
                            size_t iso_fId =
                                iso_face_id_of_tet_face[std::make_pair(tetId, face_orient.first)];
                            size_t shell_start =
                                (face_orient.second == 1)
                                    ? shell_of_half_patch[2 * patch_of_face[iso_fId]]
                                    : shell_of_half_patch[2 * patch_of_face[iso_fId] + 1];
                            // follow the ray till another iso-vertex or the sink
                            auto v_curr = extreme_v2;
                            while (next_vert[v_curr] != Arrangement<3>::None &&
                                   iso_vert_on_v_v_next[v_curr] == Arrangement<3>::None) {
                                v_curr = next_vert[v_curr];
                            }
                            if (iso_vert_on_v_v_next[v_curr] != Arrangement<3>::None) {
                                // reached iso-vert at end of the ray
                                auto iso_vId_end = iso_vert_on_v_v_next[v_curr];
                                auto end_tetId = iso_verts[iso_vId_end].tet_index;
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
                                compute_edge_intersection_order(end_tet_cut_result,
                                    local_v1,
                                    local_v2,
                                    sorted_vert_indices_on_edge);
                                // find the end shell
                                compute_passing_face(end_tet_cut_result,
                                    sorted_vert_indices_on_edge[1],
                                    sorted_vert_indices_on_edge.front(),
                                    face_orient);
                                iso_fId = iso_face_id_of_tet_face[std::make_pair(
                                    end_tetId, face_orient.first)];
                                size_t shell_end =
                                    (face_orient.second == 1)
                                        ? shell_of_half_patch[2 * patch_of_face[iso_fId]]
                                        : shell_of_half_patch[2 * patch_of_face[iso_fId] + 1];
                                // link start shell with end shell
                                shell_links.emplace_back(shell_start, shell_end);
                            } else {
                                // next_vert[v_curr] is None, v_curr is the sink vertex
                                // link shell_start with the sink
                                shell_links.emplace_back(shell_start, Arrangement<3>::None);
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
                    timing_labels.emplace_back("arrCells(group shells into arrCells)");
                    ScopedTimer<> timer("arrangement cells: group shells into arrangement cells");
                    compute_arrangement_cells(shells.size(), shell_links, arrangement_cells);
                    timings.push_back(timer.toc());
                }
                //            std::cout << "arrangement cells = " << std::endl;
                //            for (const auto& arr_cell : arrangement_cells) {
                //                std::cout << "(";
                //                for (auto s : arr_cell) {
                //                    std::cout << s << ", ";
                //                }
                //                std::cout << ")" << std::endl;
                //            }
            } else { // group simplicial cells into arrangement cells
                // ------------------- build adjacency graph of simplicial cells ---------------
                // pair (tetId, tet_cell_Id)
                std::vector<std::pair<size_t, size_t>> tet_cell_of_simp_cell;
                // simplicial half-face info
                // if the face is on isosurface, let s be its shell index, record (-s-1)
                // if the face is not on isosurface, record the simp_cell index on the opposite side
                std::vector<long long> simp_half_face_info;
                std::vector<size_t> simp_hFace_start_index;
                {
                    timing_labels.emplace_back("arrCells(build simpCell graph)");
                    ScopedTimer<> timer("arrangement cells: build simplicial cells graph");
                    size_t num_simplicial_cells = 0;
                    size_t num_simp_half_faces = 0;
                    for (size_t i = 0; i < tets.size(); i++) {
                        if (cut_result_index[i] != Arrangement<3>::None) {
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
                        if (cut_result_index[i] != Arrangement<3>::None) {
                            size_t iso_fId_start_index = iso_fId_start_index_of_tet[i];
                            size_t global_vId_start_index = global_vId_start_index_of_tet[i];
                            const auto& arrangement = cut_results[cut_result_index[i]];
                            const auto& faces = arrangement.faces;
                            const auto& cells = arrangement.cells;
                            for (size_t j = 0; j < cells.size(); j++) {
                                // create a new simplicial cell
                                const auto& cell = cells[j];
                                size_t cur_simp_cell_id = tet_cell_of_simp_cell.size();
                                tet_cell_of_simp_cell.emplace_back(i, j);
                                for (size_t k = 0; k < cell.faces.size(); k++) {
                                    size_t fId = cell.faces[k];
                                    size_t iso_fId = iso_fId_of_tet_face[iso_fId_start_index + fId];
                                    if (iso_fId !=
                                        Arrangement<3>::None) { // face k is on iso-surface
                                        size_t patch_id = patch_of_face[iso_fId];
                                        size_t half_patch_id = (faces[fId].positive_cell == j)
                                                                   ? (2 * patch_id)
                                                                   : (2 * patch_id + 1);
                                        size_t shell_id = shell_of_half_patch[half_patch_id];
                                        simp_half_face_info.push_back(-shell_id - 1);
                                    } else { // face k is not on iso-surface
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
                        } else { // tet i has no isosurface
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
                    timing_labels.emplace_back("arrCells(group simpCells into arrCells)");
                    ScopedTimer<> timer("arrangement cells: group simplicial cells");
                    size_t num_simp_cells = tet_cell_of_simp_cell.size();
                    std::vector<bool> visited_simp_cell(num_simp_cells, false);
                    std::vector<absl::flat_hash_set<size_t>> arrangement_cell_incident_shells;
                    for (size_t i = 0; i < num_simp_cells; i++) {
                        if (!visited_simp_cell[i]) {
                            // new arrangement cell
                            arrangement_cell_incident_shells.emplace_back();
                            auto& incident_shells = arrangement_cell_incident_shells.back();
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
                                        // the half-face is on isosurface
                                        size_t shell_id = -simp_half_face_info[start_index + j] - 1;
                                        incident_shells.insert(shell_id);
                                    } else {
                                        // the half-face is not on isosurface
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
                    arrangement_cells.reserve(arrangement_cell_incident_shells.size());
                    for (auto & arrangement_cell_incident_shell : arrangement_cell_incident_shells) {
                        arrangement_cells.emplace_back(arrangement_cell_incident_shell.begin(),
                            arrangement_cell_incident_shell.end());
                    }
                    timings.push_back(timer.toc());
                }
            }
        }
        timings.push_back(timer.toc());
        timing_labels.emplace_back("arrangement cells");
    }
    std::cout << "num_cells = " << arrangement_cells.size() << std::endl;
    stats_labels.emplace_back("num_cells");
    stats.push_back(arrangement_cells.size());

    if (components.size() > 1) {
        timing_labels.emplace_back("arrCells(other)");
        size_t num_timings = timings.size();
        if (use_topo_ray_shooting) {
            timings.push_back(timings[num_timings - 1] - timings[num_timings - 2] -
                              timings[num_timings - 3] - timings[num_timings - 4] -
                              timings[num_timings - 5] - timings[num_timings - 6]);
        } else {
            // baseline: group simplicial cells into arrangement cells
            timings.push_back(
                timings[num_timings - 1] - timings[num_timings - 2] - timings[num_timings - 3]);
        }
    }

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
            shells,
            components,
            arrangement_cells);
        save_result_msh(output_dir + "/iso_mesh",
            iso_pts,
            iso_faces,
            patches,
            iso_edges,
            chains,
            non_manifold_edges_of_vert,
            half_patch_list,
            shells,
            components,
            arrangement_cells);
        //
        if (components.size() > 1) {
            save_nesting_data(
                output_dir + "/nesting_data.json", next_vert, extremal_edge_of_component);
        }
    }
    // test: export timings
    save_timings(output_dir + "/timings.json", timing_labels, timings);
    // export statistics
    save_statistics(output_dir + "/stats.json", stats_labels, stats);


    return 0;
}
