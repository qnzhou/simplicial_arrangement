#include <simplicial_arrangement/simplicial_arrangement.h>
#include <simplicial_arrangement/lookup_table.h>

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

#include <absl/container/flat_hash_map.h>

using namespace simplicial_arrangement;

// a polygon from the iso-surface of implicit function
struct IsoFace
{
    // a list of polygon's vertex indices (index into some global list of vertices)
    std::vector<size_t> vert_indices;
    // the local index of this polygon in all the tets that contains it
    // each pair is (tet_Id, tet_face_Id)
    std::vector<std::pair<size_t, size_t>> tet_face_indices;
};

// vertex of isosurface: a pair of indices (tet_id, tet_vert_id)
typedef std::pair<size_t, size_t> IsoVert;

bool load_tet_mesh(const std::string& filename,
    std::vector<std::vector<double>>& pts,
    std::vector<std::vector<size_t>>& tets)
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
        pts[j].resize(3);
        for (size_t k = 0; k < 3; k++) {
            pts[j][k] = data[0][j][k].get<double>();
        }
    }
    //
    tets.resize(data[1].size());
    for (size_t j = 0; j < tets.size(); j++) {
        tets[j].resize(4);
        for (size_t k = 0; k < 4; k++) {
            tets[j][k] = data[1][j][k].get<size_t>();
        }
    }
    return true;
}

void extract_iso_mesh(const std::vector<bool>& has_isosurface,
    const std::vector<Arrangement<3>> &cut_results,
    const std::vector<std::vector<size_t>> &func_in_tet,
    const std::vector<std::vector<size_t>>& tets,
    std::vector<IsoVert>& iso_verts,
    std::vector<IsoFace>& iso_faces) {
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
                            iso_verts[num_iso_vert] = std::make_pair(i, j);
                            iso_vId_of_vert[j] = num_iso_vert;
                            ++num_iso_vert;
                        }
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
                            iso_verts[num_iso_vert] = std::make_pair(i, j);
                            iso_vId_of_vert[j] = num_iso_vert;
                            ++num_iso_vert;
                        }
                        break;
                    }
                    case 0: // in tet cell
                        iso_verts[num_iso_vert] = std::make_pair(i,j);
                        iso_vId_of_vert[j] = num_iso_vert;
                        ++num_iso_vert;
                        break;
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
                            iso_verts[num_iso_vert] = std::make_pair(i, j);
                            iso_vId_of_vert[j] = num_iso_vert;
                            ++num_iso_vert;
                        }
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

double sphere_function(const std::vector<double> &center, double r, const std::vector<double> &p) {
    return (p[0] - center[0]) * (p[0] - center[0]) + (p[1] - center[1]) * (p[1] - center[1]) +
           (p[2] - center[2]) * (p[2] - center[2]) - r * r;
}

int sign(double x) {
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
    std::vector<std::vector<double>> pts;
    std::vector<std::vector<size_t>> tets;
    load_tet_mesh(tet_mesh_file, pts, tets);
    std::cout << "tet mesh: " << pts.size() << " verts, " << tets.size() << " tets." << std::endl;

    // load implicit function values, or evaluate
    size_t n_func = 4;
    std::vector<std::vector<double>> centers(n_func);
    centers[0] = {0, 0, 0};
    centers[1] = {0.5, 0, 0};
    centers[2] = {0, 0.5, 0};
    centers[3] = {0, 0, 0.5};
    double radius = 0.5;

    std::vector<std::vector<double>> funcVals(n_func);
    for (size_t i = 0; i < n_func; i++) {
        funcVals[i].resize(pts.size());
        for (size_t j = 0; j < pts.size(); j++) {
            funcVals[i][j] = sphere_function(centers[i], radius, pts[j]);
        }
    }
    

    // function signs at vertices
    std::vector<std::vector<int>> funcSigns(n_func);
    for (size_t i = 0; i < n_func; i++) {
        funcSigns[i].resize(pts.size());
        for (size_t j = 0; j < pts.size(); j++) {
            funcSigns[i][j] = sign(funcVals[i][j]);
        }
    }

    // find functions whose iso-surfaces intersect tets
    std::vector<std::vector<size_t>> func_in_tet(tets.size());
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

    // compute arrangement in each tet (iterative plane cut)
    std::vector<bool> has_isosurface(tets.size(), false);
    std::vector<Arrangement<3>> cut_results(tets.size());
    size_t num_intersecting_tet = 0;
    for (size_t i = 0; i < tets.size(); i++) {
        const auto& func_ids = func_in_tet[i];
        if (! func_ids.empty()) {
            has_isosurface[i] = true;
            ++num_intersecting_tet;
            size_t v1 = tets[i][0];
            size_t v2 = tets[i][1];
            size_t v3 = tets[i][2];
            size_t v4 = tets[i][3];
            std::vector<Plane<double, 3>> planes(func_ids.size());
            for (size_t j = 0; j < func_ids.size(); j++) {
                size_t f_id = func_ids[j];
                planes[j] = {
                    funcVals[f_id][v1], funcVals[f_id][v2], funcVals[f_id][v3], funcVals[f_id][v4]};
            }
            //
            cut_results[i] = compute_arrangement(planes);
        }
    }
    std::cout << "num_intersecting_tet = " << num_intersecting_tet << std::endl;

    // extract arrangement mesh
    std::vector<IsoVert> iso_verts;
    std::vector<IsoFace> iso_faces;
    extract_iso_mesh(has_isosurface, cut_results, func_in_tet, tets, iso_verts, iso_faces);
    
    // test
    std::cout << "num iso-vertices = " << iso_verts.size() << std::endl;
    std::cout << "num iso-faces = " << iso_faces.size() << std::endl;
    

    return 0;
}
