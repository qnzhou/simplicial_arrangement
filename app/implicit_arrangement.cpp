#include <simplicial_arrangement/simplicial_arrangement.h>
#include <simplicial_arrangement/lookup_table.h>

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <queue>
#include <nlohmann/json.hpp>
#include "ScopedTimer.h"
#include <chrono>
#include <absl/container/flat_hash_map.h>

typedef std::chrono::duration<double> Time_duration;

using namespace simplicial_arrangement;



// a polygon from the iso-surface of implicit function
struct IsoFace
{
    // a list of polygon's vertex indices (index into some global list of iso-vertices)
    std::vector<size_t> vert_indices;
    // the local index of this polygon in all the tets that contains it
    // each pair is (tet_Id, tet_face_Id)
    std::vector<std::pair<size_t, size_t>> tet_face_indices;
    // list of indices of bounding iso-edges (index into a global list of iso-edges)
    std::vector<size_t> edge_indices;
};

// vertex of isosurface
struct IsoVert
{
    // the tet containing the IsoVert
    size_t tet_index;
    // the index of IsoVert in tet.vertices
    size_t tet_vert_index;
    // minimal simplex that contains the IsoVert
    size_t simplex_size; // 1: point, 2: edge, 3: triangle, 4: tetrahedron
    // index into a list of tet vertices
    std::array<size_t, 4> simplex_vert_indices;
    // list of implicit functions whose isosurfaces pass IsoVert (indexed into a global list of implicit functions)
    std::array<size_t, 3> func_indices;
};

struct IsoEdge
{
    size_t v1;
    size_t v2;
    // each pair is (iso_face_Id, edge_face_Id)
    // iso_face_Id: face index in the global list of iso-faces
    // edge_face_Id: edge index in the iso-face
    std::vector<std::pair<size_t, size_t>> face_edge_indices;
};

bool load_tet_mesh(const std::string& filename, 
    std::vector<std::array<double, 3>> &pts,
    std::vector<std::array<size_t, 4>> &tets)
{
    using json = nlohmann::json;
    std::ifstream fin(filename.c_str());
    if (!fin) {
        std::cout << "tet mesh file not exist!" << std::endl;
        return false;
    }
    json data;
    fin >> data;
    fin.close();
    //
    pts.resize(data[0].size());
    for (size_t j = 0; j < pts.size(); j++) {
        for (size_t k = 0; k < 3; k++) {
            pts[j][k] = data[0][j][k].get<double>();
        }
    }
    //
    tets.resize(data[1].size());
    for (size_t j = 0; j < tets.size(); j++) {
        for (size_t k = 0; k < 4; k++) {
            tets[j][k] = data[1][j][k].get<size_t>();
        }
    }
    return true;
}

bool save_result(const std::string& filename, const std::vector<std::array<double, 3>>& iso_pts,
    const std::vector<IsoFace>& iso_faces, const std::vector<std::vector<size_t>> &patches,
    const std::vector<IsoEdge>& iso_edges, const std::vector<std::vector<size_t>> &chains,
    const std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
    const std::vector<std::vector<std::pair<size_t, int>>>& half_patch_list,
    const std::vector<std::vector<size_t>>& arrangement_cells)
{
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str());
    //
    json jPts;
    for (size_t i = 0; i < iso_pts.size(); i++) {
        jPts.push_back(json(iso_pts[i]));
    }
    //
    json jFaces;
    for (size_t i = 0; i < iso_faces.size(); i++) {
        jFaces.push_back(json(iso_faces[i].vert_indices));
    }
    // 
    json jPatches;
    for (size_t i = 0; i < patches.size(); i++) {
        jPatches.push_back(json(patches[i]));
    }
    //
    json jEdges;
    for (size_t i = 0; i < iso_edges.size(); i++) {
        jEdges.push_back({iso_edges[i].v1, iso_edges[i].v2});
    }
    //
    json jChains;
    for (size_t i = 0; i < chains.size(); i++) {
        jChains.push_back(json(chains[i]));
    }
    //
    json jCorners;
    for (size_t i = 0; i < non_manifold_edges_of_vert.size(); i++) {
        if (non_manifold_edges_of_vert[i].size() > 2) {
            jCorners.push_back(i);
        }
    }
    //
    json jHalfPatchsList;
    for (size_t i = 0; i < half_patch_list.size(); i++) {
        json jHalfPatchs;
        for (size_t j = 0; j < half_patch_list[i].size(); j++) {
            jHalfPatchs.push_back(json(half_patch_list[i][j]));
        }
        jHalfPatchsList.push_back(jHalfPatchs);
    }
    //
    json jArrCells;
    for (size_t i = 0; i < arrangement_cells.size(); i++) {
        jArrCells.push_back(json(arrangement_cells[i]));
    }
    //
    json jOut;
    jOut.push_back(jPts);
    jOut.push_back(jFaces);
    jOut.push_back(jPatches);
    jOut.push_back(jEdges);
    jOut.push_back(jChains);
    jOut.push_back(jCorners);
    jOut.push_back(jHalfPatchsList);
    jOut.push_back(jArrCells);
    fout << jOut << std::endl;
    fout.close();
    return true;
}

void extract_iso_mesh(const std::vector<bool>& has_isosurface,
    const std::vector<Arrangement<3>> &cut_results,
    const std::vector<std::vector<size_t>> &func_in_tet,
    const std::vector<std::array<size_t, 4>>& tets,
    std::vector<IsoVert>& iso_verts,
    std::vector<IsoFace>& iso_faces) {
    ScopedTimer<> timer("extract iso mesh (topology only)");
    // hash table for vertices on the boundary of tetrahedron
    absl::flat_hash_map<size_t, size_t> vert_on_tetVert;
    absl::flat_hash_map<std::array<size_t, 3>, size_t> vert_on_tetEdge;
    absl::flat_hash_map<std::array<size_t, 5>, size_t> vert_on_tetFace;
    // hash table for faces on the boundary of tetrahedron
    absl::flat_hash_map<std::array<size_t, 3>, size_t> face_on_tetFace;
    //
    size_t num_tet = has_isosurface.size();
    size_t max_num_vert = 0;
    size_t max_num_face = 0;
    for (size_t i = 0; i < num_tet; i++) {
        if (has_isosurface[i]) {
            max_num_vert += cut_results[i].vertices.size();
            max_num_face += cut_results[i].faces.size();
        }
    }
    iso_verts.resize(max_num_vert);
    iso_faces.resize(max_num_face);
    //std::cout << "max_num_vert = " << max_num_vert << std::endl;
    //std::cout << "max_num_face = " << max_num_face << std::endl;
    size_t num_iso_vert = 0;
    size_t num_iso_face = 0;
    //
    for (size_t i = 0; i < num_tet; i++) {
        if (has_isosurface[i]) {
            const auto& arrangement = cut_results[i];
            const auto& vertices = arrangement.vertices;
            const auto& faces = arrangement.faces;
            const auto& func_ids = func_in_tet[i];
            auto num_func = func_ids.size();
            // find vertices and faces on isosurface
            std::vector<bool> is_iso_vert(vertices.size(), false);
            std::vector<bool> is_iso_face(faces.size(), false);
            for (size_t j = 0; j < faces.size(); j++) {
                auto pid = faces[j].supporting_plane;
                auto uid = arrangement.unique_plane_indices[pid];
                for (const auto& plane_id : arrangement.unique_planes[uid]) {
                    if (plane_id > 3) {  // plane 0,1,2,3 are tet boundaries
                        is_iso_face[j] = true;
                        for (const auto& vid : faces[j].vertices) {
                            is_iso_vert[vid] = true;
                        }
                        break;
                    }
                }
            }
            // create iso-vertices
            std::vector<size_t> iso_vId_of_vert(vertices.size(), Arrangement<3>::None);
            for (size_t j = 0; j < vertices.size(); j++) {
                if (is_iso_vert[j]) {
                    std::array<size_t, 3> implicit_pIds;
                    std::array<size_t, 3> bndry_pIds;
                    size_t num_bndry_planes = 0;
                    size_t num_impl_planes = 0;
                    const auto& vertex = vertices[j];
                    // vertex.size() == 3
                    for (size_t k = 0; k < 3; k++) {
                        if (vertex[k] > 3) { // plane 0,1,2,3 are tet boundaries
                            implicit_pIds[num_impl_planes] = func_ids[vertex[k] - 4];
                            ++num_impl_planes;
                        } else {
                            bndry_pIds[num_bndry_planes] = vertex[k];
                            ++num_bndry_planes;
                        }
                    }
                    switch (num_bndry_planes) {
                    case 2: // on tet edge
                    {
                        std::vector<bool> used_pId(4, false);
                        used_pId[bndry_pIds[0]] = true;
                        used_pId[bndry_pIds[1]] = true;
                        std::array<size_t, 2> vIds;
                        size_t num_vIds = 0;
                        for (size_t k = 0; k < 4; k++) {
                            if (!used_pId[k]) {
                                vIds[num_vIds] = tets[i][k];
                                ++num_vIds;
                            }
                        }
                        size_t vId1 = vIds[0];
                        size_t vId2 = vIds[1];
                        if (vId1 > vId2) {
                            size_t tmp = vId1;
                            vId1 = vId2;
                            vId2 = tmp;                        
                        }
                        std::array<size_t, 3> key = {vId1, vId2, implicit_pIds[0]};
                        auto iter_inserted = vert_on_tetEdge.try_emplace(key, num_iso_vert);
                        if (iter_inserted.second) {
                            auto& iso_vert = iso_verts[num_iso_vert];
                            iso_vert.tet_index = i;
                            iso_vert.tet_vert_index = j;
                            iso_vert.simplex_size = 2;
                            iso_vert.simplex_vert_indices[0] = vId1;
                            iso_vert.simplex_vert_indices[1] = vId2;
                            iso_vert.func_indices[0] = implicit_pIds[0];
                            ++num_iso_vert;
                        }
                        iso_vId_of_vert[j] = iter_inserted.first->second;
                        break;
                    }
                    case 1: // on tet face
                    {
                        std::array<size_t, 3> vIds;
                        size_t pId = bndry_pIds[0];
                        size_t num_vIds = 0;
                        for (size_t k = 0; k < 4; k++) {
                            if (k != pId) {
                                vIds[num_vIds] = tets[i][k];
                                ++num_vIds;
                            }
                        }
                        std::sort(vIds.begin(), vIds.end());
                        std::array<size_t, 5> key = {
                            vIds[0], vIds[1], vIds[2], implicit_pIds[0], implicit_pIds[1]};
                        auto iter_inserted = vert_on_tetFace.try_emplace(key, num_iso_vert);
                        if (iter_inserted.second) {
                            auto& iso_vert = iso_verts[num_iso_vert];
                            iso_vert.tet_index = i;
                            iso_vert.tet_vert_index = j;
                            iso_vert.simplex_size = 3;
                            iso_vert.simplex_vert_indices[0] = vIds[0];
                            iso_vert.simplex_vert_indices[1] = vIds[1];
                            iso_vert.simplex_vert_indices[2] = vIds[2];
                            iso_vert.func_indices[0] = implicit_pIds[0];
                            iso_vert.func_indices[1] = implicit_pIds[1];
                            ++num_iso_vert;
                        }
                        iso_vId_of_vert[j] = iter_inserted.first->second;
                        break;
                    }
                    case 0: // in tet cell
                    {
                        auto& iso_vert = iso_verts[num_iso_vert];
                        iso_vert.tet_index = i;
                        iso_vert.tet_vert_index = j;
                        iso_vert.simplex_size = 4;
                        iso_vert.simplex_vert_indices = tets[i];
                        iso_vert.func_indices = implicit_pIds;
                        iso_vId_of_vert[j] = num_iso_vert;
                        ++num_iso_vert;
                        break;
                    }
                    case 3: // on tet vertex
                    {
                        std::vector<bool> used_pId(4, false);
                        for (const auto& pId : bndry_pIds) {
                            used_pId[pId] = true;
                        }
                        size_t vId;
                        for (size_t k = 0; k < 4; k++) {
                            if (!used_pId[k]) {
                                vId = k;
                                break;
                            }
                        }
                        auto key = tets[i][vId];
                        auto iter_inserted = vert_on_tetVert.try_emplace(key, num_iso_vert);
                        if (iter_inserted.second) {
                            auto& iso_vert = iso_verts[num_iso_vert];
                            iso_vert.tet_index = i;
                            iso_vert.tet_vert_index = j;
                            iso_vert.simplex_size = 1;
                            iso_vert.simplex_vert_indices[0] = tets[i][vId];
                            ++num_iso_vert;
                        }
                        iso_vId_of_vert[j] = iter_inserted.first->second;
                        break;
                    }
                    default: break;
                    }
                }
            }            
            // create iso-faces
            for (size_t j = 0; j < faces.size(); j++) {
                if (is_iso_face[j]) {
                    std::vector<size_t> face_verts(faces[j].vertices.size());
                    for (size_t k = 0; k < face_verts.size(); k++) {
                        face_verts[k] = iso_vId_of_vert[faces[j].vertices[k]];
                    }
                    //
                    auto pid = faces[j].supporting_plane;
                    auto uid = arrangement.unique_plane_indices[pid];
                    bool face_on_tet_boundary = false;
                    for (const auto& plane_id : arrangement.unique_planes[uid]) {
                        if (plane_id < 4) { // plane 0,1,2,3 are tet boundaries
                            face_on_tet_boundary = true;
                            break;
                        }
                    }
                    if (face_on_tet_boundary) {
                        std::vector<size_t> sorted_face_verts = face_verts;
                        std::sort(sorted_face_verts.begin(), sorted_face_verts.end());
                        std::array<size_t, 3> key = {
                            sorted_face_verts[0], sorted_face_verts[1], sorted_face_verts[2]};
                        auto iter_inserted = face_on_tetFace.try_emplace(key, num_iso_face);
                        if (iter_inserted.second) {
                            iso_faces[num_iso_face].vert_indices = face_verts;
                            iso_faces[num_iso_face].tet_face_indices.emplace_back(i, j);
                            ++num_iso_face;
                        } else { // iso_face inserted before
                            size_t iso_face_id = (iter_inserted.first)->second;
                            iso_faces[iso_face_id].tet_face_indices.emplace_back(i, j);
                        }                    
                    } else { // face not on tet boundary
                        iso_faces[num_iso_face].vert_indices = face_verts;
                        iso_faces[num_iso_face].tet_face_indices.emplace_back(i, j);
                        ++num_iso_face;
                    }
                }
            }
            
        }
    }
    //
    iso_verts.resize(num_iso_vert);
    iso_faces.resize(num_iso_face);
}

// compute iso-edges and edge-face connectivity
void compute_iso_edges(std::vector<IsoFace>& iso_faces, std::vector<IsoEdge>& iso_edges) 
{
    ScopedTimer<> timer("compute iso-edges and edge-face connectivity");
    size_t max_num_edge = 0;
    for (size_t i = 0; i < iso_faces.size(); i++) {
        max_num_edge += iso_faces[i].vert_indices.size();
    }
    iso_edges.resize(max_num_edge);
    size_t num_iso_edge = 0;
    // map: (v1, v2) -> iso-edge index
    absl::flat_hash_map<std::pair<size_t, size_t>, size_t> edge_id;
    for (size_t i = 0; i < iso_faces.size(); i++) {
        auto& face = iso_faces[i];
        size_t num_edge = face.vert_indices.size();
        face.edge_indices.resize(num_edge);
        for (size_t j = 0; j < num_edge; j++) {
            size_t v1 = face.vert_indices[j];
            size_t v2 = (j+1==num_edge) ? face.vert_indices[0] : face.vert_indices[j + 1];
            // swap if v1 > v2
            size_t tmp = v1;
            if (v1 > v2) {
                v1 = v2;
                v2 = tmp;
            }
            //            
            auto iter_inserted = edge_id.try_emplace(std::make_pair(v1, v2), num_iso_edge);
            if (iter_inserted.second) { // new iso-edge
                iso_edges[num_iso_edge].v1 = v1;
                iso_edges[num_iso_edge].v2 = v2;
                iso_edges[num_iso_edge].face_edge_indices.emplace_back(i, j);
                face.edge_indices[j] = num_iso_edge; 
                ++num_iso_edge;                    
            } else { // existing iso-edge
                size_t eId = iter_inserted.first->second;
                iso_edges[eId].face_edge_indices.emplace_back(i, j);
                face.edge_indices[j] = eId;
            }            
        }
    }
    //
    iso_edges.resize(num_iso_edge);
}

// group iso-faces into patches
void compute_patches(const std::vector<IsoFace> &iso_faces, const std::vector<IsoEdge> &iso_edges,
    std::vector<std::vector<size_t>> &patches)
{
    ScopedTimer<> timer("group iso-faces into patches");
    std::vector<bool> visisted_face(iso_faces.size(), false);
    for (size_t i = 0; i < iso_faces.size(); i++) {
        if (!visisted_face[i]) {
            // new patch
            patches.emplace_back();
            auto& patch = patches.back();
            std::queue<size_t> Q;
            Q.push(i);
            patch.emplace_back(i);
            visisted_face[i] = true;
            while (! Q.empty()) {
                auto fId = Q.front();
                Q.pop();
                for (size_t eId : iso_faces[fId].edge_indices) {
                    if (iso_edges[eId].face_edge_indices.size() == 2) { // manifold edge                        
                        size_t other_fId = (iso_edges[eId].face_edge_indices[0].first == fId)
                                             ? iso_edges[eId].face_edge_indices[1].first
                                             : iso_edges[eId].face_edge_indices[0].first;
                        if (!visisted_face[other_fId]) {
                            Q.push(other_fId);
                            patch.emplace_back(other_fId);
                            visisted_face[other_fId] = true;
                        }
                    }
                }
            }
        }
    }
}

// group non-manifold iso-edges into chains
void compute_chains(const std::vector<IsoEdge>& iso_edges, 
    const std::vector<std::vector<size_t>> &non_manifold_edges_of_vert, 
    std::vector<std::vector<size_t>> &chains) 
{
    std::vector<bool> visited_edge(iso_edges.size(), false);
    for (size_t i = 0; i < iso_edges.size(); i++) {
        if (!visited_edge[i] && iso_edges[i].face_edge_indices.size() > 2) {
            // unvisited non-manifold iso-edge (not a boundary edge)
            // new chain
            chains.emplace_back();
            auto& chain = chains.back();
            std::queue<size_t> Q;
            Q.push(i);
            chain.emplace_back(i);
            visited_edge[i] = true;
            while (!Q.empty()) {
                auto eId = Q.front();
                Q.pop();
                // v1
                size_t v = iso_edges[eId].v1;
                if (non_manifold_edges_of_vert[v].size() == 2) {
                    size_t other_eId = (non_manifold_edges_of_vert[v][0] == eId)
                                           ? non_manifold_edges_of_vert[v][1]
                                           : non_manifold_edges_of_vert[v][0];
                    if (!visited_edge[other_eId]) {
                        Q.push(other_eId);
                        chain.emplace_back(other_eId);
                        visited_edge[other_eId] = true;
                    }
                }
                // v2
                v = iso_edges[eId].v2;
                if (non_manifold_edges_of_vert[v].size() == 2) {
                    size_t other_eId = (non_manifold_edges_of_vert[v][0] == eId)
                                           ? non_manifold_edges_of_vert[v][1]
                                           : non_manifold_edges_of_vert[v][0];
                    if (!visited_edge[other_eId]) {
                        Q.push(other_eId);
                        chain.emplace_back(other_eId);
                        visited_edge[other_eId] = true;
                    }
                }                
            }
        }
    }
}

// compute neighboring pair of half-faces around an iso-edge in a tetrahedron
// pair<size_t, int> : pair (iso-face index, iso-face orientation)
void compute_face_order_in_one_tet(const Arrangement<3>& tet_cut_result,
    const std::vector<IsoFace>& iso_faces,
    const IsoEdge& iso_edge,     
    std::vector<std::pair<size_t, int>>& ordered_faces)
{
    // in a single tetrahedron, num face pairs == num iso-faces incident to the iso-edge
    size_t num_incident_faces = iso_edge.face_edge_indices.size();
    ordered_faces.resize(num_incident_faces);
    //
    std::vector<bool> is_incident_faces(tet_cut_result.faces.size(), false);
    // map: tet face id --> iso face id
    std::vector<size_t> iso_face_Id_of_face(tet_cut_result.faces.size(), Arrangement<3>::None);
    for (const auto& fId_eId_pair : iso_edge.face_edge_indices) {
        size_t iso_face_id = fId_eId_pair.first;
        size_t face_id = iso_faces[iso_face_id].tet_face_indices[0].second;
        is_incident_faces[face_id] = true;
        iso_face_Id_of_face[face_id] = iso_face_id;
    }
    

    // travel around the edge
    size_t iso_face_id1 = iso_edge.face_edge_indices[0].first;
    size_t face_id1 = iso_faces[iso_face_id1].tet_face_indices[0].second;   
    int face_sign = 1;    
    size_t cell_id = tet_cut_result.faces[face_id1].positive_cell;
    size_t face_id2 = Arrangement<3>::None;
    for (size_t i = 0; i < num_incident_faces; i++) {
        ordered_faces[i] = std::make_pair(iso_face_id1, face_sign);
        // find next face        
        for (const auto& fId : tet_cut_result.cells[cell_id].faces) {
            if (is_incident_faces[fId] && fId != face_id1) {
                face_id2 = fId;
            }
        }
        // assert: face_id2 != None
        // get sign of face2 and find next cell
        if (tet_cut_result.faces[face_id2].positive_cell == cell_id) {
            cell_id = tet_cut_result.faces[face_id2].negative_cell;
            face_sign = -1;
        } else {
            cell_id = tet_cut_result.faces[face_id2].positive_cell;
            face_sign = 1;
        }
        // update face1 and clear face2
        face_id1 = face_id2;
        iso_face_id1 = iso_face_Id_of_face[face_id2];
        face_id2 = Arrangement<3>::None;
    }
}

void compute_arrangement_cells(
    size_t num_patch, const std::vector<std::vector<std::pair<size_t, int>>>& half_patch_list, 
    std::vector<std::vector<size_t>> &arrangement_cells)
{
    ScopedTimer<> timer("group patches into arrangement cells");
    // (patch i, 1) <--> 2i,  (patch i, -1) <--> 2i+1
    // compute half-patch adjacency list
    std::vector<std::vector<size_t>> half_patch_adj_list(2 * num_patch);
    for (const auto& half_patches : half_patch_list) {
        for (size_t i = 0; i+1 < half_patches.size(); i++) {
            const auto& hp1 = half_patches[i];
            const auto& hp2 = half_patches[i + 1];
            // half-patch index of hp1
            size_t hp_Id1 = (hp1.second == 1) ? 2 * hp1.first : (2 * hp1.first + 1);
            // half-patch index of hp2.opposite
            size_t hp_Id2 = (hp2.second == 1) ? (2 * hp2.first + 1) : 2 * hp2.first;
            half_patch_adj_list[hp_Id1].push_back(hp_Id2);
            half_patch_adj_list[hp_Id2].push_back(hp_Id1);
        }
        //
        const auto& hp1 = half_patches.back();
        const auto& hp2 = half_patches.front();
        size_t hp_Id1 = (hp1.second == 1) ? 2 * hp1.first : (2 * hp1.first + 1);
        size_t hp_Id2 = (hp2.second == 1) ? (2 * hp2.first + 1) : 2 * hp2.first;
        half_patch_adj_list[hp_Id1].push_back(hp_Id2);
        half_patch_adj_list[hp_Id2].push_back(hp_Id1);
    }
    // find connected component of half-patch adjacency graph
    std::vector<bool> visited_half_patch(2 * num_patch, false);
    arrangement_cells.clear();
    for (size_t i = 0; i < 2*num_patch; i++) {
        if (!visited_half_patch[i]) {
            // create new component
            arrangement_cells.emplace_back();
            auto& cell = arrangement_cells.back();
            std::queue<size_t> Q;
            Q.push(i);
            cell.push_back(i);
            visited_half_patch[i] = true;
            while (! Q.empty()) {
                auto half_patch = Q.front();
                Q.pop();
                for (auto hp : half_patch_adj_list[half_patch]) {
                    if (!visited_half_patch[hp]) {
                        cell.push_back(hp);
                        Q.push(hp);
                        visited_half_patch[hp] = true;
                    }
                }
            }
        }
    }
    // get arrangement cells as list of patch indices
    for (auto& cell : arrangement_cells) {
        for (auto& pId : cell) {
            pId /= 2;   // get patch index of half-patch
        }
    }
}

// compute barycentric coordinate of Point (intersection of three planes)
// Point in tet cell
//template <typename Scalar>
inline void compute_barycentric_coords(
    const std::array<double, 4>& plane1, const std::array<double, 4>& plane2, const std::array<double, 4>& plane3,
    std::array<double,4>& bary_coords)
{
    double n1 = plane1[3] * (plane2[2] * plane3[1] - plane2[1] * plane3[2]) +
                plane1[2] * (plane2[1] * plane3[3] - plane2[3] * plane3[1]) +
                plane1[1] * (plane2[3] * plane3[2] - plane2[2] * plane3[3]);
    double n2 = plane1[3] * (plane2[0] * plane3[2] - plane2[2] * plane3[0]) +
                plane1[2] * (plane2[3] * plane3[0] - plane2[0] * plane3[3]) +
                plane1[0] * (plane2[2] * plane3[3] - plane2[3] * plane3[2]);
    double n3 = plane1[3] * (plane2[1] * plane3[0] - plane2[0] * plane3[1]) +
                plane1[1] * (plane2[0] * plane3[3] - plane2[3] * plane3[0]) +
                plane1[0] * (plane2[3] * plane3[1] - plane2[1] * plane3[3]);
    double n4 = plane1[2] * (plane2[0] * plane3[1] - plane2[1] * plane3[0]) +
                plane1[1] * (plane2[2] * plane3[0] - plane2[0] * plane3[2]) +
                plane1[0] * (plane2[1] * plane3[2] - plane2[2] * plane3[1]);
    double d = n1 + n2 + n3 + n4;
    //
    bary_coords[0] = n1 / d;
    bary_coords[1] = n2 / d;
    bary_coords[2] = n3 / d;
    bary_coords[3] = n4 / d;
}

// Point on tet face
inline void compute_barycentric_coords(const std::array<double, 3>& plane1,
    const std::array<double, 3>& plane2,
    std::array<double, 3>& bary_coords)
{
    double n1 = plane1[2] * plane2[1] - plane1[1] * plane2[2];
    double n2 = plane1[0] * plane2[2] - plane1[2] * plane2[0];
    double n3 = plane1[1] * plane2[0] - plane1[0] * plane2[1];
    double d = n1 + n2 + n3;
    //
    bary_coords[0] = n1 / d;
    bary_coords[1] = n2 / d;
    bary_coords[2] = n3 / d;
}

// Point on tet edge
inline void compute_barycentric_coords(double f1, double f2,
    std::array<double, 2>& bary_coords)
{
    bary_coords[0] = f2 / (f2 - f1);
    bary_coords[1] = 1 - bary_coords[0];
}


// implicit functions
inline double sphere_function(const std::array<double,3> &center, double r, const std::array<double,3> &p) {
    return (p[0] - center[0]) * (p[0] - center[0]) + (p[1] - center[1]) * (p[1] - center[1]) +
           (p[2] - center[2]) * (p[2] - center[2]) - r * r;
}

inline int sign(double x) {
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

int main(int argc, const char* argv[])
{
    std::cout << "load table ..." << std::endl;
    bool loaded = load_lookup_table();
    if (loaded) {
        std::cout << "loading complete." << std::endl;
    }

    // load tet mesh
    std::string tet_mesh_file = "D:/research/simplicial_arrangement/data/tet_mesh.json";
    std::vector<std::array<double, 3>> pts;
    std::vector<std::array<size_t, 4>> tets;
    load_tet_mesh(tet_mesh_file, pts, tets);
    std::cout << "tet mesh: " << pts.size() << " verts, " << tets.size() << " tets." << std::endl;

    // load implicit function values, or evaluate
    size_t n_func = 4;
    std::vector<std::array<double,3>> centers(n_func);
    centers[0] = {0, 0, 0};
    centers[1] = {0.5, 0, 0};
    centers[2] = {0, 0.5, 0};
    centers[3] = {0, 0, 0.5};
    double radius = 0.5;

    std::vector<std::vector<double>> funcVals(n_func);
    {
        ScopedTimer<> timer("evaluate implicit functions at vertices");
        for (size_t i = 0; i < n_func; i++) {
            funcVals[i].resize(pts.size());
            for (size_t j = 0; j < pts.size(); j++) {
                funcVals[i][j] = sphere_function(centers[i], radius, pts[j]);
            }
        }
    }
    

    // function signs at vertices
    std::vector<std::vector<int>> funcSigns(n_func);
    {
        ScopedTimer<> timer("compute function signs at vertices");
        for (size_t i = 0; i < n_func; i++) {
            funcSigns[i].resize(pts.size());
            for (size_t j = 0; j < pts.size(); j++) {
                funcSigns[i][j] = sign(funcVals[i][j]);
            }
        }
    }

    // find functions whose iso-surfaces intersect tets
    std::vector<std::vector<size_t>> func_in_tet(tets.size());
    {
        ScopedTimer<> timer("filter out intersecting implicits in each tet");
        for (size_t i = 0; i < tets.size(); i++) {
            func_in_tet.reserve(n_func);
            for (size_t j = 0; j < n_func; j++) {
                int pos_count = 0;
                int neg_count = 0;
                for (const auto& vId : tets[i]) {
                    if (funcSigns[j][vId] == 1) {
                        pos_count += 1;
                    } else if (funcSigns[j][vId] == -1) {
                        neg_count += 1;
                    }
                }
                // tets[i].size() == 4
                if (pos_count < 4 && neg_count < 4) {
                    func_in_tet[i].push_back(j);
                }
            }
        }
    }

    // test: marching tet on each implicit function
    {
        ScopedTimer<> timer("test: marching tet");
        for (size_t i = 0; i < funcVals.size(); i++) {
            for (size_t j = 0; j < tets.size(); j++) {
                size_t pos_count = 0;
                size_t neg_count = 0;
                for (const auto& vId : tets[j]) {
                    if (funcSigns[i][vId] == 1) {
                        pos_count += 1;
                    } else if (funcSigns[i][vId] == -1) {
                        neg_count += 1;
                    }
                }                
                if (pos_count < 4 && neg_count < 4) {
                    // iso-surface of function i intersects tet j
                    std::vector<Plane<double, 3>> planes(1);
                    planes[0] = {funcVals[i][tets[j][0]],
                        funcVals[i][tets[j][1]],
                        funcVals[i][tets[j][2]],
                        funcVals[i][tets[j][3]]};
                    compute_arrangement(planes);
                }
                
            }
        }
    }

    // compute arrangement in each tet (iterative plane cut)
    std::vector<bool> has_isosurface(tets.size(), false);
    std::vector<Arrangement<3>> cut_results(tets.size());
    size_t num_intersecting_tet = 0;
    //
    Time_duration time_1_func = Time_duration::zero();
    Time_duration time_2_func = Time_duration::zero();
    Time_duration time_more_func = Time_duration::zero();
    //
    {
        ScopedTimer<> timer("compute arrangement in all tets");
        for (size_t i = 0; i < tets.size(); i++) {
            const auto& func_ids = func_in_tet[i];
            if (!func_ids.empty()) {
                std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
                //
                has_isosurface[i] = true;
                ++num_intersecting_tet;
                size_t v1 = tets[i][0];
                size_t v2 = tets[i][1];
                size_t v3 = tets[i][2];
                size_t v4 = tets[i][3];
                std::vector<Plane<double, 3>> planes(func_ids.size());
                for (size_t j = 0; j < func_ids.size(); j++) {
                    size_t f_id = func_ids[j];
                    planes[j] = {funcVals[f_id][v1],
                        funcVals[f_id][v2],
                        funcVals[f_id][v3],
                        funcVals[f_id][v4]};
                }
                //
                cut_results[i] = compute_arrangement(planes);
                //
                std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
                switch (func_ids.size()) {
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
    extract_iso_mesh(has_isosurface, cut_results, func_in_tet, tets, iso_verts, iso_faces);
    //std::cout << "num iso-vertices = " << iso_verts.size() << std::endl;
    //std::cout << "num iso-faces = " << iso_faces.size() << std::endl;
    
    // compute xyz of iso-vertices
    std::vector<std::array<double, 3>> iso_pts(iso_verts.size());
    {
        ScopedTimer<> timer("compute xyz of iso-vertices");
        for (size_t i = 0; i < iso_verts.size(); i++) {
            const auto& iso_vert = iso_verts[i];
            switch (iso_vert.simplex_size) {
            case 2: // on tet edge
            {
                auto vId1 = iso_vert.simplex_vert_indices[0];
                auto vId2 = iso_vert.simplex_vert_indices[1];
                auto fId = iso_vert.func_indices[0];
                auto f1 = funcVals[fId][vId1];
                auto f2 = funcVals[fId][vId2];
                std::array<double, 2> b;
                compute_barycentric_coords(f1, f2, b);
                iso_pts[i][0] = b[0] * pts[vId1][0] + b[1] * pts[vId2][0];
                iso_pts[i][1] = b[0] * pts[vId1][1] + b[1] * pts[vId2][1];
                iso_pts[i][2] = b[0] * pts[vId1][2] + b[1] * pts[vId2][2];
                break;
            }
            case 3: // on tet face
            {
                auto vId1 = iso_vert.simplex_vert_indices[0];
                auto vId2 = iso_vert.simplex_vert_indices[1];
                auto vId3 = iso_vert.simplex_vert_indices[2];
                auto fId1 = iso_vert.func_indices[0];
                auto fId2 = iso_vert.func_indices[1];
                std::array<double, 3> f1s = {
                    funcVals[fId1][vId1], funcVals[fId1][vId2], funcVals[fId1][vId3]};
                std::array<double, 3> f2s = {
                    funcVals[fId2][vId1], funcVals[fId2][vId2], funcVals[fId2][vId3]};
                std::array<double, 3> b;
                compute_barycentric_coords(f1s, f2s, b);
                iso_pts[i][0] = b[0] * pts[vId1][0] + b[1] * pts[vId2][0] + b[2] * pts[vId3][0];
                iso_pts[i][1] = b[0] * pts[vId1][1] + b[1] * pts[vId2][1] + b[2] * pts[vId3][1];
                iso_pts[i][2] = b[0] * pts[vId1][2] + b[1] * pts[vId2][2] + b[2] * pts[vId3][2];
                break;
            }
            case 4: // in tet cell
            {
                auto vId1 = iso_vert.simplex_vert_indices[0];
                auto vId2 = iso_vert.simplex_vert_indices[1];
                auto vId3 = iso_vert.simplex_vert_indices[2];
                auto vId4 = iso_vert.simplex_vert_indices[3];
                auto fId1 = iso_vert.func_indices[0];
                auto fId2 = iso_vert.func_indices[1];
                auto fId3 = iso_vert.func_indices[2];
                std::array<double, 4> f1s = {funcVals[fId1][vId1],
                    funcVals[fId1][vId2],
                    funcVals[fId1][vId3],
                    funcVals[fId1][vId4]};
                std::array<double, 4> f2s = {funcVals[fId2][vId1],
                    funcVals[fId2][vId2],
                    funcVals[fId2][vId3],
                    funcVals[fId2][vId4]};
                std::array<double, 4> f3s = {funcVals[fId3][vId1],
                    funcVals[fId3][vId2],
                    funcVals[fId3][vId3],
                    funcVals[fId3][vId4]};
                std::array<double, 4> b;
                compute_barycentric_coords(f1s, f2s, f3s, b);
                iso_pts[i][0] = b[0] * pts[vId1][0] + b[1] * pts[vId2][0] + b[2] * pts[vId3][0] +
                                b[3] * pts[vId4][0];
                iso_pts[i][1] = b[0] * pts[vId1][1] + b[1] * pts[vId2][1] + b[2] * pts[vId3][1] +
                                b[3] * pts[vId4][1];
                iso_pts[i][2] = b[0] * pts[vId1][2] + b[1] * pts[vId2][2] + b[2] * pts[vId3][2] +
                                b[3] * pts[vId4][2];
                break;
            }
            case 1: // on tet vertex
                iso_pts[i] = pts[iso_vert.simplex_vert_indices[0]];
                break;
            default: break;
            }
        }
    }
    
    //  compute iso-edges and edge-face connectivity
    std::vector<IsoEdge> iso_edges;
    compute_iso_edges(iso_faces, iso_edges);
    //std::cout << "num iso-edges = " << iso_edges.size() << std::endl;

    // group iso-faces into patches
    std::vector<std::vector<size_t>> patches;
    compute_patches(iso_faces, iso_edges, patches);
    //std::cout << "num patches = " << patches.size() << std::endl;    

    
    // group non-manifold iso-edges into chains
    std::vector<std::vector<size_t>> non_manifold_edges_of_vert(iso_pts.size());
    std::vector<std::vector<size_t>> chains;
    {
        ScopedTimer<> timer("group non-manifold iso-edges into chains");
        // get incident non-manifold edges for iso-vertices
        for (size_t i = 0; i < iso_edges.size(); i++) {
            if (iso_edges[i].face_edge_indices.size() > 2) 
            { // non-manifold edge (not a boundary edge)
                // there is only one patch indicent to a boundary edge,
                // so there is no need to figure out the "order" of patches around a boundary edge
                non_manifold_edges_of_vert[iso_edges[i].v1].push_back(i);
                non_manifold_edges_of_vert[iso_edges[i].v2].push_back(i);
            }
        }
        // group non-manifold iso-edges into chains
        compute_chains(iso_edges, non_manifold_edges_of_vert, chains);
        // std::cout << "num chains = " << chains.size() << std::endl;
    }

    // compute list of incident tets for each vertex
    std::vector<std::vector<size_t>> incident_tets_of_vert(pts.size());
    {
        ScopedTimer<> timer("compute vertex-tet adjacency of input tet mesh");
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
    std::vector<std::vector<std::pair<size_t, int>>> half_faces_list(chains.size());
    std::vector<std::vector<std::pair<size_t, int>>> half_patch_list(half_faces_list.size());
    {
        ScopedTimer<> timer("compute order of patches around chains");
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

        // compute map: iso-face Id --> patch Id
        std::vector<size_t> patch_of_face(iso_faces.size());
        for (size_t i = 0; i < patches.size(); i++) {
            for (const auto& fId : patches[i]) {
                patch_of_face[fId] = i;
            }
        }

        // replace iso-face indices by patch indices
        // std::vector<std::vector<std::pair<size_t, int>>> half_patch_list(half_faces_list.size());
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
    std::vector<std::vector<size_t>> arrangement_cells;
    compute_arrangement_cells(patches.size(), half_patch_list, arrangement_cells);
    //std::cout << "num arrangement cells = " << arrangement_cells.size() << std::endl;

    // test: export iso-mesh, patches, chains
    save_result("D:/research/simplicial_arrangement/data/iso_mesh.json",
        iso_pts,
        iso_faces,
        patches,
        iso_edges,
        chains,
        non_manifold_edges_of_vert,
        half_patch_list,
        arrangement_cells);
    
    return 0;
}
