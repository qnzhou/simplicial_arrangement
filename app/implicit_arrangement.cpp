#include <simplicial_arrangement/simplicial_arrangement.h>
#include <simplicial_arrangement/lookup_table.h>

#include <iostream>
#include <fstream>
#include <queue>
#include <nlohmann/json.hpp>
#include "ScopedTimer.h"
#include <chrono>
#include <absl/container/flat_hash_map.h>
#include <absl/container/flat_hash_set.h>

#include <Eigen/Core>


#include "implicit_arrangement_util.h"

typedef std::chrono::duration<double> Time_duration;

using namespace simplicial_arrangement;



int main(int argc, const char* argv[])
{
    // choose method to obtain arrangement cell
    // 1. order patches around chains
    // 2. group simplicial cells into arrangement cells
    bool use_group_simplicial_cells_into_arrangement_cells = false;

    std::cout << "load table ..." << std::endl;
    bool loaded = load_lookup_table();
    if (loaded) {
        std::cout << "loading complete." << std::endl;
    }

    // load tet mesh
    //    std::string dataDir = "D:/research/simplicial_arrangement/data/";
    std::string dataDir = "/Users/charlesdu/Downloads/research/implicit_modeling/code/simplicial_arrangement/data/";
    std::string resolution = "100k";
    std::string tet_mesh_file = dataDir + "tet_mesh_" + resolution + ".json";
    std::vector<std::array<double, 3>> pts;
    std::vector<std::array<size_t, 4>> tets;
    load_tet_mesh(tet_mesh_file, pts, tets);
    size_t n_tets = tets.size();
    size_t n_pts = pts.size();
    std::cout << "tet mesh: " << pts.size() << " verts, " << tets.size() << " tets." << std::endl;

    // find all boundary triangles of the tet mesh
    {
        ScopedTimer<> timer("find boundary triangles of tet mesh");
        absl::flat_hash_set<std::array<size_t, 3>> unpaired_tris;
        std::vector<std::array<size_t, 3>> four_tris(4);
        for (size_t i = 0; i < tets.size(); i++) {
            auto tet = tets[i];
            std::sort(tet.begin(), tet.end());
            four_tris[0] = {tet[1], tet[2], tet[3]};
            four_tris[1] = {tet[0], tet[2], tet[3]};
            four_tris[2] = {tet[0], tet[1], tet[3]};
            four_tris[3] = {tet[0], tet[1], tet[2]};
            for (size_t j = 0; j < 4; j++) {
                auto iter_inserted = unpaired_tris.insert(four_tris[j]);
                if (!iter_inserted.second) {
                    // triangle inserted before, we found a pair of triangles, delete it
                    unpaired_tris.erase(iter_inserted.first);
                }
            }
        }
        std::cout << "num boundary tris = " << unpaired_tris.size() << std::endl;
    }

    // find all boundary triangles of the tet mesh
//    {
//        ScopedTimer<> timer("find boundary triangles of tet mesh (std::move)");
//        absl::flat_hash_set<std::array<size_t, 3>> unpaired_tris;
//        std::vector<std::array<size_t, 3>> four_tris(4);
//        for (size_t i = 0; i < tets.size(); i++) {
//            auto tet = tets[i];
//            std::sort(tet.begin(), tet.end());
//            four_tris[0] = {tet[1], tet[2], tet[3]};
//            four_tris[1] = {tet[0], tet[2], tet[3]};
//            four_tris[2] = {tet[0], tet[1], tet[3]};
//            four_tris[3] = {tet[0], tet[1], tet[2]};
//            for (size_t j = 0; j < 4; j++) {
//                auto iter_inserted = unpaired_tris.insert(std::move(four_tris[j]));
//                if (!iter_inserted.second) {
//                    // triangle inserted before, we found a pair of triangles, delete it
//                    unpaired_tris.erase(iter_inserted.first);
//                }
//            }
//        }
//        std::cout << "num boundary tris = " << unpaired_tris.size() << std::endl;
//    }

    // find all boundary triangles of the tet mesh
//    {
//        ScopedTimer<> timer("find boundary triangles of tet mesh (move + reserve)");
//        absl::flat_hash_set<std::array<size_t, 3>> unpaired_tris;
//        unpaired_tris.reserve(n_tets * 4);
//        std::vector<std::array<size_t, 3>> four_tris(4);
//        for (size_t i = 0; i < tets.size(); i++) {
//            auto tet = tets[i];
//            std::sort(tet.begin(), tet.end());
//            four_tris[0] = {tet[1], tet[2], tet[3]};
//            four_tris[1] = {tet[0], tet[2], tet[3]};
//            four_tris[2] = {tet[0], tet[1], tet[3]};
//            four_tris[3] = {tet[0], tet[1], tet[2]};
//            for (size_t j = 0; j < 4; j++) {
//                auto iter_inserted = unpaired_tris.insert(std::move(four_tris[j]));
//                if (!iter_inserted.second) {
//                    // triangle inserted before, we found a pair of triangles, delete it
//                    unpaired_tris.erase(iter_inserted.first);
//                }
//            }
//        }
//        std::cout << "num boundary tris = " << unpaired_tris.size() << std::endl;
//    }




    // extract tet boundary mesh
    //    std::vector<std::array<double, 3>> boundary_pts;
    //    std::vector<std::array<size_t, 3>> boundary_tris;
    //    extract_tet_boundary_mesh(pts, tets, boundary_pts, boundary_tris);
    //    save_tri_mesh(dataDir + "tet_bndry_mesh_" + resolution + ".json",
    //        boundary_pts, boundary_tris);


    // load implicit function values, or evaluate
    size_t n_func = 4;
    std::vector<std::array<double,3>> centers(n_func);
    centers[0] = {0, 0, 0};
    centers[1] = {0.5, 0, 0};
    centers[2] = {0, 0.5, 0};
    centers[3] = {0, 0, 0.5};
    double radius = 0.5;


    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> funcVals;
    {
        ScopedTimer<> timer("evaluate implicit functions at vertices");
        funcVals.resize(n_func, n_pts);
        for (size_t i = 0; i < n_func; i++) {
            for (size_t j = 0; j < n_pts; j++) {
                funcVals(i,j) = sphere_function(centers[i], radius, pts[j]);
            }
        }
    }


    // function signs at vertices
    Eigen::MatrixXi funcSigns;
    {
        ScopedTimer<> timer("compute function signs at vertices");
        funcSigns.resize(n_func, n_pts);
        for (size_t i = 0; i < n_func; i++) {
            for (size_t j = 0; j < n_pts; j++) {
                funcSigns(i,j) = sign(funcVals(i,j));
            }
        }
    }
//    {
//        ScopedTimer<> timer("compute function signs at vertices");
//        funcSigns = funcVals.unaryExpr([](double x) { return sign(x); });
//    }


    // find functions whose iso-surfaces intersect tets
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> func_in_tet;
    Eigen::VectorXi num_func_in_tet;
    {
        ScopedTimer<> timer("filter out intersecting implicits in each tet");
        func_in_tet.resize(n_tets, n_func);
        num_func_in_tet.setZero(n_tets);
        int pos_count;
        int neg_count;
        int num_func;
        for (size_t i = 0; i < n_tets; i++) {
            num_func = 0;
            for (size_t j = 0; j < n_func; j++) {
                pos_count = 0;
                neg_count = 0;
                for (const auto& vId : tets[i]) {
                    if (funcSigns(j,vId) == 1) {
                        pos_count += 1;
                    } else if (funcSigns(j,vId) == -1) {
                        neg_count += 1;
                    }
                }
                // tets[i].size() == 4
                if (pos_count < 4 && neg_count < 4) {
                    func_in_tet(i, num_func) = j;
                    ++num_func;
                }
            }
            num_func_in_tet[i] = num_func;
        }
    }

    // compute arrangement in each tet (iterative plane cut)
    std::vector<bool> has_isosurface;
    std::vector<Arrangement<3>> cut_results;
    size_t num_intersecting_tet = 0;
    //
    Time_duration time_1_func = Time_duration::zero();
    Time_duration time_2_func = Time_duration::zero();
    Time_duration time_more_func = Time_duration::zero();
    //
    {
        ScopedTimer<> timer("compute arrangement in all tets");
        has_isosurface.resize(n_tets, false);
        cut_results.resize(n_tets);
        int num_func;
        for (size_t i = 0; i < tets.size(); i++) {
            const auto& func_ids = func_in_tet.row(i);
            num_func = num_func_in_tet(i);
            if (num_func != 0) {
                std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
                //
                has_isosurface[i] = true;
                ++num_intersecting_tet;
                size_t v1 = tets[i][0];
                size_t v2 = tets[i][1];
                size_t v3 = tets[i][2];
                size_t v4 = tets[i][3];
                std::vector<Plane<double, 3>> planes(num_func);
                for (size_t j = 0; j < num_func; j++) {
                    size_t f_id = func_ids(j);
                    planes[j] = {funcVals(f_id,v1),
                        funcVals(f_id,v2),
                        funcVals(f_id,v3),
                        funcVals(f_id,v4)};
                }
                //
                if (num_func == 2) {
                    disable_lookup_table();
                    cut_results[i] = compute_arrangement(planes);
                    enable_lookup_table();
                } else {
                    cut_results[i] = compute_arrangement(planes);
                }

                //
                std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
                switch (num_func) {
                case 1: time_1_func += std::chrono::duration_cast<Time_duration>(t2 - t1); break;
                case 2: time_2_func += std::chrono::duration_cast<Time_duration>(t2 - t1); break;
                default: time_more_func += std::chrono::duration_cast<Time_duration>(t2 - t1); break;
                }
            }
        }
    }
    //std::cout << "num_intersecting_tet = " << num_intersecting_tet << std::endl;
    std::cout << "-- [1   func]: " << time_1_func.count() << "s" << std::endl;
    std::cout << "-- [2   func]: " << time_2_func.count() << "s" << std::endl;
    std::cout << "-- [>=3 func]: " << time_more_func.count() << "s" << std::endl;


    // extract arrangement mesh
    std::vector<IsoVert> iso_verts;
    std::vector<IsoFace> iso_faces;
    std::vector<std::vector<size_t>> global_vId_of_tet_vert;  // this is used later when we want to know the global vert index of a local vertex in a tet
    std::vector<std::vector<size_t>> iso_fId_of_tet_face; // this is used later when we want to know the global face index of a local face in a tet
    if (use_group_simplicial_cells_into_arrangement_cells) {
        extract_iso_mesh(has_isosurface,
            cut_results,
            func_in_tet,
            num_func_in_tet,
            tets,
            iso_verts,
            iso_faces,
            global_vId_of_tet_vert,
            iso_fId_of_tet_face);
    } else {
        extract_iso_mesh_pure(has_isosurface, cut_results,
            func_in_tet, num_func_in_tet, tets, iso_verts, iso_faces);
    }
    
    //std::cout << "num iso-vertices = " << iso_verts.size() << std::endl;
    //std::cout << "num iso-faces = " << iso_faces.size() << std::endl;
    
    // compute xyz of iso-vertices
    std::vector<std::array<double, 3>> iso_pts;
    compute_iso_vert_xyz(iso_verts, funcVals, pts, iso_pts);
    
    //  compute iso-edges and edge-face connectivity
    std::vector<IsoEdge> iso_edges;
    compute_iso_edges(iso_faces, iso_edges);
    //std::cout << "num iso-edges = " << iso_edges.size() << std::endl;

    // group iso-faces into patches
    std::vector<std::vector<size_t>> patches;
    compute_patches(iso_faces, iso_edges, patches);
    //std::cout << "num patches = " << patches.size() << std::endl;    

    // compute map: iso-face Id --> patch Id
    std::vector<size_t> patch_of_face;
    {
        ScopedTimer<> timer("compute iso-face-id to patch-id map");
        patch_of_face.resize(iso_faces.size());
        for (size_t i = 0; i < patches.size(); i++) {
            for (const auto& fId : patches[i]) {
                patch_of_face[fId] = i;
            }
        }
    }

    std::vector<std::vector<size_t>> arrangement_cells;
    if (use_group_simplicial_cells_into_arrangement_cells)
    // baseline approach: group simplicial cells into arrangement cells
    {
        {
            ScopedTimer<> timer("assign global index to all vertices");
            // ------------ assign a global index for every vertex (iso-vertices or tet vertices)
            // ---------------- global index of iso-vertex = iso-vertex's index in 'iso_verts'
            // global index of tet vertex = iso_verts.size() + vertex's index in 'pts'
            /*size_t tmp_count_1 = 0;
            size_t tmp_count_2 = 0;
            size_t tmp_count_3 = 0;*/
            size_t num_iso_vert = iso_verts.size();
            for (size_t i = 0; i < tets.size(); i++) {
                if (has_isosurface[i]) {
                    auto& global_vId_of_vert = global_vId_of_tet_vert[i];
                    for (size_t j = 0; j < global_vId_of_vert.size(); j++) {
                        if (global_vId_of_vert[j] == Arrangement<3>::None) {
                            global_vId_of_vert[j] = num_iso_vert + tets[i][j];
                            /*if (j > 4) {
                                ++tmp_count_1;
                            } else {
                                const auto &plane_Ids = cut_results[i].vertices[j];
                                std::vector<bool> used_pId(4, false);
                                for (const auto& pId : plane_Ids) {
                                    used_pId[pId] = true;
                                }
                                size_t vId;
                                for (size_t k = 0; k < 4; k++) {
                                    if (!used_pId[k]) {
                                        vId = k;
                                        break;
                                    }
                                }
                                if (vId == j) {
                                    ++tmp_count_2;
                                } else {
                                    ++tmp_count_3;
                                }
                            }*/
                        }
                    }

                } else { // tet i has no iso-surface
                    global_vId_of_tet_vert[i][0] = num_iso_vert + tets[i][0];
                    global_vId_of_tet_vert[i][1] = num_iso_vert + tets[i][1];
                    global_vId_of_tet_vert[i][2] = num_iso_vert + tets[i][2];
                    global_vId_of_tet_vert[i][3] = num_iso_vert + tets[i][3];
                }
            }
            // test
            // assert: tmp_count_1 == 0, tmp_count_3 == 0
            /*std::cout << "count(j >  4) = " << tmp_count_1 << std::endl;
            std::cout << "count(j <= 4 and j==vId) = " << tmp_count_2 << std::endl;
            std::cout << "count(j <= 4 and j!=vId) = " << tmp_count_3 << std::endl;*/
        }
        // ------------------- build adjacency graph of simplicial cells ---------------
        std::vector<Simplicial_Cell> simp_cells;
        {
            ScopedTimer<> timer("build adjacency graph of simplicial cells");
            size_t num_simplicial_cells = 0;
            for (size_t i = 0; i < tets.size(); i++) {
                if (has_isosurface[i]) {
                    num_simplicial_cells += cut_results[i].cells.size();
                } else {
                    ++num_simplicial_cells;
                }
            }
            simp_cells.resize(num_simplicial_cells);
            // hash table: face key --> (simplicial cell index, face index in simplicial cell)
            absl::flat_hash_map<std::array<size_t, 3>, std::pair<size_t, size_t>>
                incident_cell_of_face;
            size_t cur_simp_cell_id = 0;
            for (size_t i = 0; i < tets.size(); i++) {
                if (has_isosurface[i]) {
                    const auto& arrangement = cut_results[i];
                    const auto& cells = arrangement.cells;
                    for (size_t j = 0; j < cells.size(); j++) {
                        // create a new simplicial cell
                        const auto& cell = cells[j];
                        auto& cur_simp_cell = simp_cells[cur_simp_cell_id];
                        cur_simp_cell.tet_Id = i;
                        cur_simp_cell.tet_cell_Id = j;
                        cur_simp_cell.is_iso_face.resize(cell.faces.size(), false);
                        cur_simp_cell.face_info.resize(cell.faces.size(), Arrangement<3>::None);
                        for (size_t k = 0; k < cell.faces.size(); k++) {
                            size_t fId = cell.faces[k];
                            size_t iso_fId = iso_fId_of_tet_face[i][fId];
                            if (iso_fId != Arrangement<3>::None) { // face k is on iso-surface
                                cur_simp_cell.is_iso_face[k] = true;
                                cur_simp_cell.face_info[k] = patch_of_face[iso_fId];
                            } else { // face k is not on iso-surface
                                std::vector<size_t> face_verts = arrangement.faces[fId].vertices;
                                // convert from local vert index to global vert index
                                for (size_t l = 0; l < face_verts.size(); l++) {
                                    face_verts[l] = global_vId_of_tet_vert[i][face_verts[l]];
                                }
                                //
                                std::array<size_t, 3> key;
                                compute_iso_face_key(face_verts, key);
                                auto iter_inserted = incident_cell_of_face.try_emplace(
                                    key, std::make_pair(cur_simp_cell_id, k));
                                if (!iter_inserted.second) { // same face has been inserted before
                                    size_t opposite_simp_cell_id =
                                        iter_inserted.first->second.first;
                                    size_t opposite_cell_face_id =
                                        iter_inserted.first->second.second;
                                    // make neighbor
                                    cur_simp_cell.face_info[k] = opposite_simp_cell_id;
                                    simp_cells[opposite_simp_cell_id]
                                        .face_info[opposite_cell_face_id] = cur_simp_cell_id;
                                    // delete face in hash table since a face can only be shared by two cells 
                                    incident_cell_of_face.erase(iter_inserted.first); // this actually makes it slower
                                }
                            }
                        }
                        ++cur_simp_cell_id;
                    }
                } else { // tet i has no isosurface
                    auto& cur_simp_cell = simp_cells[cur_simp_cell_id];
                    cur_simp_cell.tet_Id = i;
                    cur_simp_cell.tet_cell_Id = 0;
                    cur_simp_cell.is_iso_face.resize(
                        4, false); // 4 triangles of tet i are not on isosurface
                    cur_simp_cell.face_info.resize(4, Arrangement<3>::None);
                    //
                    const auto& global_vIds = global_vId_of_tet_vert[i];
                    std::vector<std::array<size_t, 3>> four_face_verts(4);
                    // face 0
                    four_face_verts[0] = {global_vIds[1], global_vIds[2], global_vIds[3]};
                    std::sort(four_face_verts[0].begin(), four_face_verts[0].end());
                    // face 1
                    four_face_verts[1] = {global_vIds[2], global_vIds[3], global_vIds[0]};
                    std::sort(four_face_verts[1].begin(), four_face_verts[1].end());
                    // face 2
                    four_face_verts[2] = {global_vIds[3], global_vIds[0], global_vIds[1]};
                    std::sort(four_face_verts[2].begin(), four_face_verts[2].end());
                    // face 3
                    four_face_verts[3] = {global_vIds[0], global_vIds[1], global_vIds[2]};
                    std::sort(four_face_verts[3].begin(), four_face_verts[3].end());
                    //
                    for (size_t j = 0; j < 4; j++) {
                        auto iter_inserted = incident_cell_of_face.try_emplace(
                            four_face_verts[j], std::make_pair(cur_simp_cell_id, j));
                        if (!iter_inserted.second) { // same face has been inserted before
                            size_t opposite_simp_cell_id = iter_inserted.first->second.first;
                            size_t opposite_cell_face_id = iter_inserted.first->second.second;
                            // make neighbor
                            cur_simp_cell.face_info[j] = opposite_simp_cell_id;
                            simp_cells[opposite_simp_cell_id].face_info[opposite_cell_face_id] =
                                cur_simp_cell_id;
                            // delete face in hash table since a face can only be shared by two cells
                            incident_cell_of_face.erase(iter_inserted.first); // this actually makes it slower
                        }
                    }
                    ++cur_simp_cell_id;
                }
            }
            //std::cout << "incident_cell_of_face.size() = " << incident_cell_of_face.size() << std::endl;
        }
        {
            // ------------------- group simplical cells into arrangement cells ---------------
            ScopedTimer<> timer("group simplical cells into arrangement cells");
            std::vector<bool> visited_simp_cell(simp_cells.size(), false);
            std::vector<std::vector<size_t>> simp_cells_of_arrangement_cell;
            // arrangement_cell_incident_patch[i][j] == true if cell i is incident to patch j
            std::vector<std::vector<bool>> arrangement_cell_incident_patch;
            for (size_t i = 0; i < simp_cells.size(); i++) {
                if (!visited_simp_cell[i]) {
                    // new arrangement cell
                    simp_cells_of_arrangement_cell.emplace_back();
                    auto& simp_cell_list = simp_cells_of_arrangement_cell.back();
                    arrangement_cell_incident_patch.emplace_back(patches.size(), false);
                    auto& incident_patches = arrangement_cell_incident_patch.back();
                    //
                    std::queue<size_t> Q;
                    Q.push(i);
                    visited_simp_cell[i] = true;
                    while (!Q.empty()) {
                        size_t simp_cell_id = Q.front();
                        //simp_cell_list.push_back(simp_cell_id);  // note: maybe this is not needed
                        Q.pop();
                        const auto& simp_cell = simp_cells[simp_cell_id];
                        for (size_t j = 0; j < simp_cell.face_info.size(); j++) {
                            if (simp_cell.is_iso_face[j]) {
                                incident_patches[simp_cell.face_info[j]] = true;
                            } else if (simp_cell.face_info[j] != Arrangement<3>::None &&
                                       !visited_simp_cell[simp_cell.face_info[j]]) {
                                // if the opposite cell exists and has not been visited
                                size_t opposite_simp_cell_id = simp_cell.face_info[j];
                                Q.push(opposite_simp_cell_id);
                                visited_simp_cell[opposite_simp_cell_id] = true;
                            }
                        }
                    }
                }
            }
            // collect bounding patches for each arrangement cell
            arrangement_cells.resize(arrangement_cell_incident_patch.size());
            for (size_t i = 0; i < arrangement_cells.size(); i++) {
                auto& arrangement_cell = arrangement_cells[i];
                const auto& incident_patches = arrangement_cell_incident_patch[i];
                for (size_t j = 0; j < incident_patches.size(); j++) {
                    if (incident_patches[j]) {
                        arrangement_cell.push_back(j);
                    }
                }
            }
            // std::cout << "num arrangement cells = " << arrangement_cells.size() << std::endl;
        }
    }
    // another approach: order patches around chains
    else {
        // group non-manifold iso-edges into chains
        std::vector<std::vector<size_t>> non_manifold_edges_of_vert;
        std::vector<std::vector<size_t>> chains;
        {
            ScopedTimer<> timer("group non-manifold iso-edges into chains");
            non_manifold_edges_of_vert.resize(iso_pts.size());
            // get incident non-manifold edges for iso-vertices
            for (size_t i = 0; i < iso_edges.size(); i++) {
                if (iso_edges[i].face_edge_indices.size() >
                    2) { // non-manifold edge (not a boundary edge)
                    // there is only one patch indicent to a boundary edge,
                    // so there is no need to figure out the "order" of patches around a boundary
                    // edge
                    non_manifold_edges_of_vert[iso_edges[i].v1].push_back(i);
                    non_manifold_edges_of_vert[iso_edges[i].v2].push_back(i);
                }
            }
            // group non-manifold iso-edges into chains
            compute_chains(iso_edges, non_manifold_edges_of_vert, chains);
            // std::cout << "num chains = " << chains.size() << std::endl;
        }

        // compute list of incident tets for each vertex
        std::vector<std::vector<size_t>> incident_tets_of_vert;
        {
            ScopedTimer<> timer("compute vertex-tet adjacency of input tet mesh");
            incident_tets_of_vert.resize(pts.size());
            // compute number of tets incident to each vertex
            std::vector<size_t> num_incident_tets(pts.size(), 0);
            for (size_t i = 0; i < tets.size(); i++) {
                const auto& tet = tets[i];
                num_incident_tets[tet[0]] += 1;
                num_incident_tets[tet[1]] += 1;
                num_incident_tets[tet[2]] += 1;
                num_incident_tets[tet[3]] += 1;
            }
            // compute indices of tets incident to each vertex
            std::vector<size_t> size_incident_tets(pts.size(), 0);
            for (size_t i = 0; i < pts.size(); i++) {
                incident_tets_of_vert[i].resize(num_incident_tets[i]);
            }
            for (size_t i = 0; i < tets.size(); i++) {
                const auto& tet = tets[i];
                incident_tets_of_vert[tet[0]][size_incident_tets[tet[0]]] = i;
                incident_tets_of_vert[tet[1]][size_incident_tets[tet[1]]] = i;
                incident_tets_of_vert[tet[2]][size_incident_tets[tet[2]]] = i;
                incident_tets_of_vert[tet[3]][size_incident_tets[tet[3]]] = i;
                size_incident_tets[tet[0]] += 1;
                size_incident_tets[tet[1]] += 1;
                size_incident_tets[tet[2]] += 1;
                size_incident_tets[tet[3]] += 1;
            }
        }

        // compute order of patches around chains
        std::vector<std::vector<std::pair<size_t, int>>> half_faces_list;
        std::vector<std::vector<std::pair<size_t, int>>> half_patch_list;
        {
            ScopedTimer<> timer("compute order of patches around chains");
            half_faces_list.resize(chains.size());
            half_patch_list.resize(half_faces_list.size());
            // pick representative iso-edge from each chain
            std::vector<size_t> chain_representatives(chains.size());
            for (size_t i = 0; i < chains.size(); i++) {
                chain_representatives[i] = chains[i][0];
            }
            // order iso-faces incident to each representative iso-edge
            // pair<size_t, int> : pair (iso-face index, iso-face orientation)
            /*std::vector<std::vector<std::pair<size_t, int>>> half_faces_list(chains.size());*/
            for (size_t i = 0; i < chain_representatives.size(); i++) {
                // first try: assume each representative edge is in the interior of a tetrahedron
                const auto& iso_edge = iso_edges[chain_representatives[i]];
                auto iso_face_id = iso_edge.face_edge_indices[0].first;
                auto tet_id = iso_faces[iso_face_id].tet_face_indices[0].first;
                compute_face_order_in_one_tet(
                    cut_results[tet_id], iso_faces, iso_edge, half_faces_list[i]);
            }


            // replace iso-face indices by patch indices
            // std::vector<std::vector<std::pair<size_t, int>>>
            // half_patch_list(half_faces_list.size());
            for (size_t i = 0; i < half_faces_list.size(); i++) {
                half_patch_list[i].resize(half_faces_list[i].size());
                for (size_t j = 0; j < half_faces_list[i].size(); j++) {
                    half_patch_list[i][j] = std::make_pair(
                        patch_of_face[half_faces_list[i][j].first], half_faces_list[i][j].second);
                }
            }
        }

        // group patches into arrangement cells
        // each cell is a represented as a list of patch indices
        //std::vector<std::vector<size_t>> arrangement_cells;
        compute_arrangement_cells(patches.size(), half_patch_list, arrangement_cells);
        std::cout << "num arrangement cells = " << arrangement_cells.size() << std::endl;

        // test: export iso-mesh, patches, chains
        save_result(
            dataDir + "iso_mesh_" + resolution + ".json",
            iso_pts,
            iso_faces,
            patches,
            iso_edges,
            chains,
            non_manifold_edges_of_vert,
            half_patch_list,
            arrangement_cells);
    }

    if (use_group_simplicial_cells_into_arrangement_cells) {
        save_result_mini(
            dataDir + "iso_mesh_" + resolution + ".json",
            iso_pts,
            iso_faces,
            patches,
            arrangement_cells);
    }
    
    
    return 0;
}
