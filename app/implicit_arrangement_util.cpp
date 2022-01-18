#include "implicit_arrangement_util.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <queue>

#include "ScopedTimer.h"
#include "io.h"

#include <absl/container/flat_hash_map.h>
#include <absl/container/flat_hash_set.h>


bool parse_config_file(const std::string& filename,
    std::string& tet_mesh_file,
    std::string& sphere_file,
    std::string& output_dir,
    bool& use_lookup,
    bool& use_2func_lookup,
    bool& use_bbox,
    std::array<double,3> &bbox_min,
    std::array<double,3> &bbox_max)
{
    using json = nlohmann::json;
    std::ifstream fin(filename.c_str());
    if (!fin) {
        std::cout << "configure file not exist!" << std::endl;
        return false;
    }
    json data;
    fin >> data;
    fin.close();
    //
    tet_mesh_file = data["tetMeshFile"];
    sphere_file = data["sphereFile"];
    output_dir = data["outputDir"];
    use_lookup = data["useLookup"];
    use_2func_lookup = data["use2funcLookup"];
    if (data.contains("useBBox")) {
        use_bbox = data["useBBox"];
    }
    if (data.contains("bboxMin")) {
        bbox_min = data["bboxMin"];
    }
    if (data.contains("bboxMax")) {
        bbox_max = data["bboxMax"];
    }
    return true;
}


bool load_tet_mesh(const std::string& filename,
    std::vector<std::array<double, 3>>& pts,
    std::vector<std::array<size_t, 4>>& tets)
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

bool load_spheres(const std::string& filename, std::vector<Sphere>& spheres)
{
    using json = nlohmann::json;
    std::ifstream fin(filename.c_str());
    if (!fin) {
        std::cout << "spheres file not exist!" << std::endl;
        return false;
    }
    json data;
    fin >> data;
    fin.close();
    //
    spheres.resize(data.size());
    for (size_t i = 0; i < spheres.size(); i++) {
        // center: {x,y,z}
        for (int j = 0; j < 3; ++j) {
            spheres[i].first[j] = data[i][0][j].get<double>();
        }
        // radius
        spheres[i].second = data[i][1];
    }
    return true;
}

bool save_result(const std::string& filename,
    const std::vector<std::array<double, 3>>& iso_pts,
    const std::vector<PolygonFace>& iso_faces,
    const std::vector<std::vector<size_t>>& patches,
    const std::vector<Edge>& edges,
    const std::vector<std::vector<size_t>>& chains,
    const std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
    const std::vector<std::vector<std::pair<size_t, int>>>& half_patch_list,
    const std::vector<std::vector<size_t>>& shells,
    const std::vector<std::vector<size_t>>& components,
    const std::vector<std::vector<size_t>>& arrangement_cells)
{
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str());
    //
    json jPts;
    for (const auto& iso_pt : iso_pts) {
        jPts.push_back(json(iso_pt));
    }
    //
    json jFaces;
    for (const auto& iso_face : iso_faces) {
        jFaces.push_back(json(iso_face.vert_indices));
    }
    //
    json jPatches;
    for (const auto& patch : patches) {
        jPatches.push_back(json(patch));
    }
    //
    json jEdges;
    for (const auto& edge : edges) {
        jEdges.push_back({edge.v1, edge.v2});
    }
    //
    json jChains;
    for (const auto& chain : chains) {
        jChains.push_back(json(chain));
    }
    //
    json jCorners;
    for (size_t i = 0; i < non_manifold_edges_of_vert.size(); i++) {
        if (non_manifold_edges_of_vert[i].size() > 2 || non_manifold_edges_of_vert[i].size() == 1) {
            jCorners.push_back(i);
        }
    }
    //
    json jHalfPatchsList;
    for (const auto& i : half_patch_list) {
        json jHalfPatchs;
        for (const auto& j : i) {
            jHalfPatchs.push_back(json(j));
        }
        jHalfPatchsList.push_back(jHalfPatchs);
    }
    //
    json jShells;
    for (const auto & shell : shells) {
        jShells.push_back(json(shell));
    }
    //
    json jComponents;
    for (const auto & component : components) {
        jComponents.push_back(json(component));
    }
    //
    json jArrCells;
    for (const auto& arrangement_cell : arrangement_cells) {
        jArrCells.push_back(json(arrangement_cell));
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
    jOut.push_back(jShells);
    jOut.push_back(jComponents);
    jOut.push_back(jArrCells);
    fout << jOut << std::endl;
    fout.close();
    return true;
}

bool save_result_msh(const std::string& filename,
    const std::vector<std::array<double, 3>>& iso_pts,
    const std::vector<PolygonFace>& iso_faces,
    const std::vector<std::vector<size_t>>& patches,
    const std::vector<Edge>& iso_edges,
    const std::vector<std::vector<size_t>>& chains,
    const std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
    const std::vector<std::vector<std::pair<size_t, int>>>& half_patch_list,
    const std::vector<std::vector<size_t>>& shells,
    const std::vector<std::vector<size_t>>& components,
    const std::vector<std::vector<size_t>>& arrangement_cells)
{
    constexpr size_t INVALID = std::numeric_limits<size_t>::max();
    std::vector<size_t> vertex_map(iso_pts.size(), INVALID);
    wmtk::MshData msh, msh2;

    // Save chains.
    auto extract_chain = [&](size_t i) {
        std::fill(vertex_map.begin(), vertex_map.end(), INVALID);
        const auto& chain = chains[i];
        size_t num_vertices = 0;
        for (auto eid : chain) {
            const auto& e = iso_edges[eid];
            if (vertex_map[e.v1] == INVALID) {
                vertex_map[e.v1] = num_vertices;
                num_vertices++;
            }
            if (vertex_map[e.v2] == INVALID) {
                vertex_map[e.v2] = num_vertices;
                num_vertices++;
            }
        }

        std::vector<std::array<double, 3>> vertices(num_vertices);
        for (size_t j=0; j<vertex_map.size(); j++) {
            if (vertex_map[j]  == INVALID) continue;
            vertices[vertex_map[j]] = iso_pts[j];
        }

        msh.add_edge_vertices(num_vertices, [&](size_t j) { return vertices[j]; });
        msh.add_edges(chain.size(), [&](size_t j) -> std::array<size_t, 2> {
                const auto eid = chain[j];
                const auto& e = iso_edges[eid];
                return {vertex_map[e.v1], vertex_map[e.v2]};
                });
    };
    for (size_t i = 0; i < chains.size(); i++) {
        extract_chain(i);
    }

    // Save patches
    auto extract_patch = [&](size_t i) -> std::tuple<std::vector<std::array<double, 3>>,
                                           std::vector<std::array<size_t, 3>>,
                                           std::vector<size_t>> {
        std::fill(vertex_map.begin(), vertex_map.end(), INVALID);
        std::vector<std::array<double, 3>> vertices;
        std::vector<std::array<size_t, 3>> triangles;
        std::vector<size_t> polygon_ids;

        size_t num_vertices = 0;
        auto add_vertex = [&](size_t vi) {
            if (vertex_map[vi] == INVALID) {
                vertex_map[vi] = num_vertices;
                num_vertices++;
            }
        };

        const auto& patch = patches[i];
        triangles.reserve(patch.size() * 2);
        polygon_ids.reserve(patch.size() * 2);

        for (const auto poly_id : patch) {
            const auto& f = iso_faces[poly_id];
            const size_t s = f.vert_indices.size();
            add_vertex(f.vert_indices[0]);
            add_vertex(f.vert_indices[1]);
            for (size_t i = 2; i < s; i++) {
                triangles.push_back({f.vert_indices[0], f.vert_indices[i - 1], f.vert_indices[i]});
                polygon_ids.push_back(poly_id);
                add_vertex(f.vert_indices[i]);
            }
        }

        vertices.resize(num_vertices);
        for (size_t j = 0; j < iso_pts.size(); j++) {
            if (vertex_map[j] == INVALID) continue;
            const auto& p = iso_pts[j];
            vertices[vertex_map[j]] = p;
        }

        for (auto& t : triangles) {
            t[0] = vertex_map[t[0]];
            t[1] = vertex_map[t[1]];
            t[2] = vertex_map[t[2]];
        }
        return {vertices, triangles, polygon_ids};
    };

    std::vector<size_t> patch_ids;
    std::vector<size_t> polygon_ids;
    for (size_t i=0; i<patches.size(); i++) {
        auto r = extract_patch(i);
        const auto& vertices = std::get<0>(r);
        const auto& triangles = std::get<1>(r);
        const auto& local_polygon_ids = std::get<2>(r);
        msh.add_face_vertices(vertices.size(), [&](size_t j) { return vertices[j]; });
        msh.add_faces(
            triangles.size(), [&](size_t j) { return triangles[j]; });
        patch_ids.insert(patch_ids.end(), triangles.size(), i);
        polygon_ids.insert(polygon_ids.end(), local_polygon_ids.begin(), local_polygon_ids.end());
    }

    msh.add_face_attribute<1>("patch_id", [&](size_t i) { return patch_ids[i]; });
    msh.add_face_attribute<1>("polygon_id", [&](size_t i) { return polygon_ids[i]; });
    msh.save(filename + "_patches.msh");

    auto extract_cell = [&](size_t i) {
        const auto& cell = arrangement_cells[i];

        std::fill(vertex_map.begin(), vertex_map.end(), INVALID);
        std::vector<std::array<double, 3>> vertices;
        std::vector<std::array<size_t, 3>> triangles;

        size_t num_vertices = 0;
        auto add_vertex = [&](size_t vi) {
            if (vertex_map[vi] == INVALID) {
                vertex_map[vi] = num_vertices;
                num_vertices++;
            }
        };

        for (auto shell_id: cell) {
            const auto& shell = shells[shell_id];
            for (auto half_patch_id : shell) {
                size_t patch_id = half_patch_id / 2;
                const auto& patch = patches[patch_id];
                for (auto poly_id : patch) {
                    const auto& f = iso_faces[poly_id];
                    const size_t s = f.vert_indices.size();
                    add_vertex(f.vert_indices[0]);
                    add_vertex(f.vert_indices[1]);
                    for (size_t i = 2; i < s; i++) {
                        triangles.push_back({f.vert_indices[0], f.vert_indices[i - 1], f.vert_indices[i]});
                        add_vertex(f.vert_indices[i]);
                    }
                }
            }
        }

        vertices.resize(num_vertices);
        for (size_t j = 0; j < iso_pts.size(); j++) {
            if (vertex_map[j] == INVALID) continue;
            const auto& p = iso_pts[j];
            vertices[vertex_map[j]] = p;
        }

        for (auto& t : triangles) {
            t[0] = vertex_map[t[0]];
            t[1] = vertex_map[t[1]];
            t[2] = vertex_map[t[2]];
        }

        msh2.add_face_vertices(num_vertices, [&](size_t j) { return vertices[j]; });
        msh2.add_faces(triangles.size(), [&](size_t j) { return triangles[j]; });
        return triangles.size();
    };

    std::vector<size_t> cell_ids;
    for (size_t i=0; i< arrangement_cells.size(); i++) {
        size_t num_faces = extract_cell(i);
        cell_ids.insert(cell_ids.end(), num_faces, i);
    }
    msh2.add_face_attribute<1>("cell_id", [&](size_t i) { return cell_ids[i]; });
    msh2.save(filename + "_cells.msh");
    return true;
}

bool save_nesting_data(const std::string& filename,
    const std::vector<size_t>& next_vert,
    const std::vector<size_t>& extremal_edge_of_component)
{
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str());
    //
    json jNextVert(next_vert);
    //
    json jExtremeEdge;
    size_t num_components = extremal_edge_of_component.size()/2;
    for (size_t i = 0; i < num_components; ++i) {
        jExtremeEdge.push_back({extremal_edge_of_component[2*i],
            extremal_edge_of_component[2*i + 1]});
    }
    //
    json jOut;
    jOut.push_back(jNextVert);
    jOut.push_back(jExtremeEdge);
    fout << jOut << std::endl;
    fout.close();
    return true;
}

bool save_result_mini(const std::string& filename,
    const std::vector<std::array<double, 3>>& iso_pts,
    const std::vector<PolygonFace>& iso_faces,
    const std::vector<std::vector<size_t>>& patches,
    const std::vector<std::vector<size_t>>& arrangement_cells)
{
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str());
    //
    json jPts;
    for (const auto& iso_pt : iso_pts) {
        jPts.push_back(json(iso_pt));
    }
    //
    json jFaces;
    for (const auto& iso_face : iso_faces) {
        jFaces.push_back(json(iso_face.vert_indices));
    }
    //
    json jPatches;
    for (const auto& patch : patches) {
        jPatches.push_back(json(patch));
    }
    //
    json jArrCells;
    for (const auto& arrangement_cell : arrangement_cells) {
        jArrCells.push_back(json(arrangement_cell));
    }
    //
    json jOut;
    jOut.push_back(jPts);
    jOut.push_back(jFaces);
    jOut.push_back(jPatches);
    jOut.push_back(jArrCells);
    fout << jOut << std::endl;
    fout.close();
    return true;
}

bool save_iso_mesh_list(const std::string& filename,
    const std::vector<std::vector<std::array<double, 3>>>& iso_pts_list,
    const std::vector<std::vector<PolygonFace>>& iso_faces_list)
{
    // assert iso_pts_list.size() == iso_faces_list.size()
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str());
    //
    json jOut;
    for (size_t i = 0; i < iso_pts_list.size(); i++) {
        const auto& iso_pts = iso_pts_list[i];
        const auto& iso_faces = iso_faces_list[i];
        json jPts;
        for (const auto& iso_pt : iso_pts) {
            jPts.push_back(json(iso_pt));
        }
        //
        json jFaces;
        for (const auto& iso_face : iso_faces) {
            jFaces.push_back(json(iso_face.vert_indices));
        }
        //
        json jMesh = {jPts, jFaces};
        //
        jOut.push_back(jMesh);
    }
    fout << jOut << std::endl;
    fout.close();
    return true;
}

bool save_tri_mesh(const std::string& filename,
    const std::vector<std::array<double, 3>>& verts,
    const std::vector<std::array<size_t, 3>>& tris)
{
    // assert iso_pts_list.size() == iso_faces_list.size()
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str());
    //
    json jVerts;
    for (const auto& vert : verts) {
        jVerts.push_back(json(vert));
    }
    //
    json jTris;
    for (const auto& tri : tris) {
        jTris.push_back(json(tri));
    }
    //
    json jMesh = {jVerts, jTris};
    fout << jMesh << std::endl;
    fout.close();
    return true;
}

bool save_timings(const std::string& filename,
    const std::vector<std::string>& timing_labels,
    const std::vector<double>& timings)
{
    // assert timing_labels.size() == timings.size()
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str());
    //
    json jOut;
    for (size_t i = 0; i < timings.size(); ++i) {
        jOut[timing_labels[i]] = timings[i];
    }
    //
    fout << jOut << std::endl;
    fout.close();
    return true;
}




// extract the boundary triangle mesh of a tet mesh
// assume: the tet mesh represents a simply-connected 3D volume
// assume: the triangles of the boundary mesh don't need to be consistently oriented
// input:
//  verts: num_vert * 3, tet mesh vertices
//  tets: num_tet * 4, tet's corner vertices' indices
// output:
//  boundary_verts: num_boundary_vert * 3, boundary vertices
//  boundary_faces: num_boundary_face * 3, boundary triangle's vertex indices
void extract_tet_boundary_mesh(const std::vector<std::array<double, 3>>& verts,
    const std::vector<std::array<size_t, 4>>& tets,
    std::vector<std::array<double, 3>>& boundary_verts,
    std::vector<std::array<size_t, 3>>& boundary_faces)
{
    ScopedTimer<> timer("extract boundary triangle mesh of tet mesh");
    absl::flat_hash_set<std::array<size_t, 3>> unpaired_tris;
    std::vector<std::array<size_t, 3>> four_tris(4);
    for (auto tet : tets) {
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
    // std::cout << "num boundary tris = " << unpaired_tris.size() << std::endl;
    boundary_faces.clear();
    boundary_faces.insert(boundary_faces.end(), unpaired_tris.begin(), unpaired_tris.end());
    // map: tet vert index --> boundary vert index
    std::vector<size_t> bndry_vId_of_tet_vert(verts.size(), Arrangement<3>::None);
    size_t num_bndry_vert = 0;
    for (auto& tri : boundary_faces) {
        for (size_t i = 0; i < 3; ++i) {
            size_t vId = tri[i];
            if (bndry_vId_of_tet_vert[vId] == Arrangement<3>::None) {
                // create new boundary vertex
                bndry_vId_of_tet_vert[vId] = num_bndry_vert;
                ++num_bndry_vert;
            }
            tri[i] = bndry_vId_of_tet_vert[vId];
        }
    }
    // copy boundary vertices
    boundary_verts.resize(num_bndry_vert);
    for (size_t i = 0; i < verts.size(); ++i) {
        if (bndry_vId_of_tet_vert[i] != Arrangement<3>::None) {
            boundary_verts[bndry_vId_of_tet_vert[i]] = verts[i];
        }
    }
}


bool save_tri_mesh_list(const std::string& filename,
    const std::vector<std::vector<std::array<double, 3>>>& verts_list,
    const std::vector<std::vector<std::array<size_t, 3>>>& tris_list)
{
    // assert verts_list.size() == tris_list.size()
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str());
    //
    json jOut;
    for (size_t i = 0; i < verts_list.size(); i++) {
        const auto& verts = verts_list[i];
        const auto& tris = tris_list[i];
        json jVerts;
        for (const auto& vert : verts) {
            jVerts.push_back(json(vert));
        }
        //
        json jTris;
        for (const auto& tri : tris) {
            jTris.push_back(json(tri));
        }
        //
        json jMesh = {jVerts, jTris};
        //
        jOut.push_back(jMesh);
    }
    fout << jOut << std::endl;
    fout.close();
    return true;
}

// given the list of vertex indices of a face, return the unique key of the face: (the smallest vert Id,
// second-smallest vert Id, the largest vert Id) assume: face_verts is a list of non-duplicate natural
// numbers, with at least three elements.
template <typename IndexType>
void compute_iso_face_key(const std::vector<IndexType>& face_verts, std::array<IndexType, 3>& key)
{
    IndexType min_vert = face_verts[0];
    size_t min_pos = 0;
    IndexType max_vert = face_verts[0];
    for (size_t i = 1; i < face_verts.size(); i++) {
        if (face_verts[i] < min_vert) {
            min_vert = face_verts[i];
            min_pos = i;
        } else if (face_verts[i] > max_vert) {
            max_vert = face_verts[i];
        }
    }
    IndexType second_min_vert = max_vert + 1;
    for (size_t i = 0; i < face_verts.size(); i++) {
        if (i != min_pos && face_verts[i] < second_min_vert) {
            second_min_vert = face_verts[i];
        }
    }
    //
    key[0] = min_vert;
    key[1] = second_min_vert;
    key[2] = max_vert;
}

void extract_iso_mesh(const std::vector<bool>& has_isosurface,
    const std::vector<Arrangement<3>>& cut_results,
    const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& func_in_tet,
    const Eigen::VectorXi& num_func_in_tet,
    const std::vector<std::array<size_t, 4>>& tets,
    std::vector<IsoVert>& iso_verts,
    std::vector<PolygonFace>& iso_faces,
    std::vector<std::vector<size_t>>& global_vId_of_tet_vert,
    std::vector<std::vector<size_t>>& iso_fId_of_tet_face)
{
    ScopedTimer<> timer("extract mesh");
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
    // std::cout << "max_num_vert = " << max_num_vert << std::endl;
    // std::cout << "max_num_face = " << max_num_face << std::endl;
    size_t num_iso_vert = 0;
    size_t num_iso_face = 0;
    // initialize map: local vert index --> iso-vert index
    global_vId_of_tet_vert.resize(num_tet);
    for (size_t i = 0; i < num_tet; i++) {
        if (has_isosurface[i]) {
            global_vId_of_tet_vert[i].resize(cut_results[i].vertices.size(), Arrangement<3>::None);
        } else {
            global_vId_of_tet_vert[i].resize(4, Arrangement<3>::None);
        }
    }
    // initialize map: local face index --> iso-face index
    iso_fId_of_tet_face.resize(num_tet);
    for (size_t i = 0; i < num_tet; i++) {
        if (has_isosurface[i]) {
            iso_fId_of_tet_face[i].resize(cut_results[i].faces.size(), Arrangement<3>::None);
        }
    }
    //
    for (size_t i = 0; i < num_tet; i++) {
        if (has_isosurface[i]) {
            const auto& arrangement = cut_results[i];
            const auto& vertices = arrangement.vertices;
            const auto& faces = arrangement.faces;
            const auto& func_ids = func_in_tet.row(i);
            auto num_func = num_func_in_tet(i);
            // find vertices and faces on isosurface
            std::vector<bool> is_iso_vert(vertices.size(), false);
            std::vector<bool> is_iso_face(faces.size(), false);
            if (arrangement.unique_planes.empty()) { // all planes are unique
                for (size_t j = 0; j < faces.size(); ++j) {
                    if (faces[j].supporting_plane > 3) { // plane 0,1,2,3 are tet boundaries
                        is_iso_face[j] = true;
                        for (const auto& vid : faces[j].vertices) {
                            is_iso_vert[vid] = true;
                        }
                    }
                }
            } else {
                for (size_t j = 0; j < faces.size(); j++) {
                    auto pid = faces[j].supporting_plane;
                    auto uid = arrangement.unique_plane_indices[pid];
                    for (const auto& plane_id : arrangement.unique_planes[uid]) {
                        if (plane_id > 3) { // plane 0,1,2,3 are tet boundaries
                            is_iso_face[j] = true;
                            for (const auto& vid : faces[j].vertices) {
                                is_iso_vert[vid] = true;
                            }
                            break;
                        }
                    }
                }
            }
            // map: local vert index --> iso-vert index
            auto& iso_vId_of_vert = global_vId_of_tet_vert[i];
            // std::vector<size_t> iso_vId_of_vert(vertices.size(), Arrangement<3>::None);
            //  create iso-vertices
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
                            implicit_pIds[num_impl_planes] = func_ids(vertex[k] - 4);
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
                    // face is on tet boundary if face.negative_cell is NONE
                    bool face_on_tet_boundary = (faces[j].negative_cell == Arrangement<3>::None);
                    //
                    if (face_on_tet_boundary) {
                        /*std::vector<size_t> sorted_face_verts = face_verts;
                        std::sort(sorted_face_verts.begin(), sorted_face_verts.end());
                        std::array<size_t, 3> key = {
                            sorted_face_verts[0], sorted_face_verts[1], sorted_face_verts[2]};*/
                        std::array<size_t, 3> key;
                        compute_iso_face_key(face_verts, key);
                        auto iter_inserted = face_on_tetFace.try_emplace(key, num_iso_face);
                        if (iter_inserted.second) {
                            iso_faces[num_iso_face].vert_indices = face_verts;
                            iso_faces[num_iso_face].tet_face_indices.emplace_back(i, j);
                            iso_fId_of_tet_face[i][j] = num_iso_face;
                            ++num_iso_face;
                        } else { // iso_face inserted before
                            size_t iso_face_id = (iter_inserted.first)->second;
                            iso_faces[iso_face_id].tet_face_indices.emplace_back(i, j);
                            iso_fId_of_tet_face[i][j] = iso_face_id;
                        }
                    } else { // face not on tet boundary
                        iso_faces[num_iso_face].vert_indices = face_verts;
                        iso_faces[num_iso_face].tet_face_indices.emplace_back(i, j);
                        iso_fId_of_tet_face[i][j] = num_iso_face;
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

void extract_iso_mesh_pure(size_t num_1_func,
    size_t num_2_func,
    size_t num_more_func,
    const std::vector<Arrangement<3>>& cut_results,
    const std::vector<size_t>& cut_result_index,
    const std::vector<size_t>& func_in_tet,
    const std::vector<size_t>& start_index_of_tet,
    const std::vector<std::array<size_t, 4>>& tets,
    std::vector<IsoVert>& iso_verts,
    std::vector<PolygonFace>& iso_faces)
{
    size_t n_tets = tets.size();
    // estimate number of iso-verts and iso-faces
    size_t max_num_face = num_1_func + 4 * num_2_func + 8 * num_more_func;
    size_t max_num_vert = max_num_face;
    iso_verts.reserve(max_num_vert);
    iso_faces.reserve(max_num_face);
    // hash table for vertices on the boundary of tetrahedron
    absl::flat_hash_map<size_t, size_t> vert_on_tetVert;
    absl::flat_hash_map<std::array<size_t, 3>, size_t> vert_on_tetEdge;
    vert_on_tetEdge.reserve(num_1_func + 3 * num_2_func + num_more_func);
    absl::flat_hash_map<std::array<size_t, 5>, size_t> vert_on_tetFace;
    vert_on_tetFace.reserve(num_2_func + 6 * num_more_func);
    // hash table for faces on the boundary of tetrahedron
    absl::flat_hash_map<std::array<size_t, 3>, size_t> face_on_tetFace;
    //
    std::vector<bool> is_iso_vert;
    is_iso_vert.reserve(8);
    std::vector<bool> is_iso_face;
    is_iso_face.reserve(9);
    std::vector<size_t> iso_vId_of_vert;
    iso_vId_of_vert.reserve(8);
    std::vector<size_t> face_verts;
    face_verts.reserve(4);
    std::array<size_t, 3> key3;
    std::array<size_t, 5> key5;
    std::array<bool, 4> used_pId;
    std::array<size_t, 2> vIds2;
    std::array<size_t, 3> vIds3;
    std::array<size_t, 3> implicit_pIds;
    std::array<size_t, 3> bndry_pIds;
    //
    for (size_t i = 0; i < n_tets; i++) {
        if (cut_result_index[i] != Arrangement<3>::None) {
            const auto& arrangement = cut_results[cut_result_index[i]];
            const auto& vertices = arrangement.vertices;
            const auto& faces = arrangement.faces;
            auto start_index = start_index_of_tet[i];
            auto num_func = start_index_of_tet[i + 1] - start_index;
            // find vertices and faces on isosurface
            is_iso_vert.clear();
            for (int j = 0; j < vertices.size(); ++j) {
                is_iso_vert.push_back(false);
            }
            is_iso_face.clear();
            if (arrangement.unique_planes.empty()) { // all planes are unique
                for (size_t j = 0; j < faces.size(); ++j) {
                    is_iso_face.push_back(false);
                    if (faces[j].supporting_plane > 3) { // plane 0,1,2,3 are tet boundaries
                        is_iso_face.back() = true;
                        for (const auto& vid : faces[j].vertices) {
                            is_iso_vert[vid] = true;
                        }
                    }
                }
            } else {
                for (size_t j = 0; j < faces.size(); j++) {
                    is_iso_face.push_back(false);
                    auto pid = faces[j].supporting_plane;
                    auto uid = arrangement.unique_plane_indices[pid];
                    for (const auto& plane_id : arrangement.unique_planes[uid]) {
                        if (plane_id > 3) { // plane 0,1,2,3 are tet boundaries
                            //                        is_iso_face[j] = true;
                            is_iso_face.back() = true;
                            for (const auto& vid : faces[j].vertices) {
                                is_iso_vert[vid] = true;
                            }
                            break;
                        }
                    }
                }
            }
            // map: local vert index --> iso-vert index
            //            std::vector<size_t> iso_vId_of_vert(vertices.size(),
            //            Arrangement<3>::None);
            iso_vId_of_vert.clear();
            // create iso-vertices
            for (size_t j = 0; j < vertices.size(); j++) {
                iso_vId_of_vert.push_back(Arrangement<3>::None);
                if (is_iso_vert[j]) {
                    //                    std::array<size_t, 3> implicit_pIds;
                    //                    std::array<size_t, 3> bndry_pIds;
                    size_t num_bndry_planes = 0;
                    size_t num_impl_planes = 0;
                    const auto& vertex = vertices[j];
                    // vertex.size() == 3
                    for (size_t k = 0; k < 3; k++) {
                        if (vertex[k] > 3) { // plane 0,1,2,3 are tet boundaries
                            //                            implicit_pIds[num_impl_planes] =
                            //                            func_ids(vertex[k] - 4);
                            implicit_pIds[num_impl_planes] =
                                func_in_tet[vertex[k] - 4 + start_index];
                            ++num_impl_planes;
                        } else {
                            bndry_pIds[num_bndry_planes] = vertex[k];
                            ++num_bndry_planes;
                        }
                    }
                    switch (num_bndry_planes) {
                    case 2: // on tet edge
                    {
                        //                        std::vector<bool> used_pId(4, false);
                        //                        std::array<bool, 4> used_pId {false, false, false,
                        //                        false}; used_pId = {false, false, false, false};
                        used_pId[0] = false;
                        used_pId[1] = false;
                        used_pId[2] = false;
                        used_pId[3] = false;
                        used_pId[bndry_pIds[0]] = true;
                        used_pId[bndry_pIds[1]] = true;
                        //                        std::array<size_t, 2> vIds;
                        size_t num_vIds = 0;
                        for (size_t k = 0; k < 4; k++) {
                            if (!used_pId[k]) {
                                vIds2[num_vIds] = tets[i][k];
                                ++num_vIds;
                            }
                        }
                        size_t vId1 = vIds2[0];
                        size_t vId2 = vIds2[1];
                        if (vId1 > vId2) {
                            size_t tmp = vId1;
                            vId1 = vId2;
                            vId2 = tmp;
                        }
                        //                        std::array<size_t, 3> key = {vId1, vId2,
                        //                        implicit_pIds[0]}; key3 = {vId1, vId2,
                        //                        implicit_pIds[0]};
                        key3[0] = vId1;
                        key3[1] = vId2;
                        key3[2] = implicit_pIds[0];
                        auto iter_inserted = vert_on_tetEdge.try_emplace(key3, iso_verts.size());
                        if (iter_inserted.second) {
                            iso_verts.emplace_back();
                            auto& iso_vert = iso_verts.back();
                            iso_vert.tet_index = i;
                            iso_vert.tet_vert_index = j;
                            iso_vert.simplex_size = 2;
                            iso_vert.simplex_vert_indices[0] = vId1;
                            iso_vert.simplex_vert_indices[1] = vId2;
                            iso_vert.func_indices[0] = implicit_pIds[0];
                        }
                        //                        iso_vId_of_vert[j] = iter_inserted.first->second;
                        iso_vId_of_vert.back() = iter_inserted.first->second;
                        break;
                    }
                    case 1: // on tet face
                    {
                        //                        std::array<size_t, 3> vIds;
                        size_t pId = bndry_pIds[0];
                        size_t num_vIds = 0;
                        for (size_t k = 0; k < 4; k++) {
                            if (k != pId) {
                                vIds3[num_vIds] = tets[i][k];
                                ++num_vIds;
                            }
                        }
                        std::sort(vIds3.begin(), vIds3.end());
                        //                        std::array<size_t, 5> key = {
                        //                            vIds[0], vIds[1], vIds[2], implicit_pIds[0],
                        //                            implicit_pIds[1]};
                        //                        key5 = {
                        //                            vIds[0], vIds[1], vIds[2], implicit_pIds[0],
                        //                            implicit_pIds[1]};
                        key5[0] = vIds3[0];
                        key5[1] = vIds3[1];
                        key5[2] = vIds3[2];
                        key5[3] = implicit_pIds[0];
                        key5[4] = implicit_pIds[1];
                        auto iter_inserted = vert_on_tetFace.try_emplace(key5, iso_verts.size());
                        if (iter_inserted.second) {
                            iso_verts.emplace_back();
                            auto& iso_vert = iso_verts.back();
                            iso_vert.tet_index = i;
                            iso_vert.tet_vert_index = j;
                            iso_vert.simplex_size = 3;
                            iso_vert.simplex_vert_indices[0] = vIds3[0];
                            iso_vert.simplex_vert_indices[1] = vIds3[1];
                            iso_vert.simplex_vert_indices[2] = vIds3[2];
                            iso_vert.func_indices[0] = implicit_pIds[0];
                            iso_vert.func_indices[1] = implicit_pIds[1];
                        }
                        //                        iso_vId_of_vert[j] = iter_inserted.first->second;
                        iso_vId_of_vert.back() = iter_inserted.first->second;
                        break;
                    }
                    case 0: // in tet cell
                    {
                        //                        iso_vId_of_vert[j] = iso_verts.size();
                        iso_vId_of_vert.back() = iso_verts.size();
                        iso_verts.emplace_back();
                        auto& iso_vert = iso_verts.back();
                        iso_vert.tet_index = i;
                        iso_vert.tet_vert_index = j;
                        iso_vert.simplex_size = 4;
                        iso_vert.simplex_vert_indices = tets[i];
                        iso_vert.func_indices = implicit_pIds;
                        break;
                    }
                    case 3: // on tet vertex
                    {
                        //                        std::vector<bool> used_pId(4, false);
                        //                        used_pId = {false, false, false, false};
                        used_pId[0] = false;
                        used_pId[1] = false;
                        used_pId[2] = false;
                        used_pId[3] = false;
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
                        auto iter_inserted = vert_on_tetVert.try_emplace(key, iso_verts.size());
                        if (iter_inserted.second) {
                            iso_verts.emplace_back();
                            auto& iso_vert = iso_verts.back();
                            iso_vert.tet_index = i;
                            iso_vert.tet_vert_index = j;
                            iso_vert.simplex_size = 1;
                            iso_vert.simplex_vert_indices[0] = tets[i][vId];
                        }
                        //                        iso_vId_of_vert[j] = iter_inserted.first->second;
                        iso_vId_of_vert.back() = iter_inserted.first->second;
                        break;
                    }
                    default: break;
                    }
                }
            }
            // create iso-faces
            for (size_t j = 0; j < faces.size(); j++) {
                if (is_iso_face[j]) {
                    //                    std::vector<size_t> face_verts(faces[j].vertices.size());
                    face_verts.clear();
                    for (size_t k = 0; k < faces[j].vertices.size(); k++) {
                        //                        face_verts[k] =
                        //                        iso_vId_of_vert[faces[j].vertices[k]];
                        face_verts.push_back(iso_vId_of_vert[faces[j].vertices[k]]);
                    }
                    //
                    // face is on tet boundary if face.negative_cell is NONE
                    bool face_on_tet_boundary = (faces[j].negative_cell == Arrangement<3>::None);
                    //
                    if (face_on_tet_boundary) {
                        /*std::vector<size_t> sorted_face_verts = face_verts;
                        std::sort(sorted_face_verts.begin(), sorted_face_verts.end());
                        std::array<size_t, 3> key = {
                            sorted_face_verts[0], sorted_face_verts[1], sorted_face_verts[2]};*/
                        //                        std::array<size_t, 3> key;
                        compute_iso_face_key(face_verts, key3);
                        auto iter_inserted = face_on_tetFace.try_emplace(key3, iso_faces.size());
                        if (iter_inserted.second) {
                            iso_faces.emplace_back();
                            iso_faces.back().vert_indices = face_verts;
                            iso_faces.back().tet_face_indices.emplace_back(i, j);
                        } else { // iso_face inserted before
                            size_t iso_face_id = (iter_inserted.first)->second;
                            iso_faces[iso_face_id].tet_face_indices.emplace_back(i, j);
                        }
                    } else { // face not on tet boundary
                        iso_faces.emplace_back();
                        iso_faces.back().vert_indices = face_verts;
                        iso_faces.back().tet_face_indices.emplace_back(i, j);
                    }
                }
            }
        }
    }
    //
}

void extract_iso_mesh_marching_tet(const std::vector<bool>& has_isosurface,
    const std::vector<Arrangement<3>>& cut_results,
    const std::vector<std::array<size_t, 4>>& tets,
    std::vector<IsoVert>& iso_verts,
    std::vector<PolygonFace>& iso_faces)
{
    // hash table for vertices on the boundary of tetrahedron
    // tet vertex index --> iso-vertex index
    absl::flat_hash_map<size_t, size_t> vert_on_tetVert;
    // (tet vertex1 index, tet vertex2 index) --> iso-vertex index
    absl::flat_hash_map<std::array<size_t, 2>, size_t> vert_on_tetEdge;
    // hash table for faces on the boundary of tetrahedron
    // (smallest iso-vert Id, second smallest iso-vert Id, largest iso-vert Id) --> iso-face index
    absl::flat_hash_map<std::array<size_t, 3>, size_t> face_on_tetFace;
    //
    size_t num_tet = has_isosurface.size();
    size_t num_non_empty_tet = 0;
    for (size_t i = 0; i < num_tet; i++) {
        if (has_isosurface[i]) {
            ++num_non_empty_tet;
        }
    }
    size_t max_num_vert = 4 * num_non_empty_tet;
    size_t max_num_face = num_non_empty_tet;
    iso_verts.resize(max_num_vert);
    iso_faces.resize(max_num_face);
    // std::cout << "max_num_vert = " << max_num_vert << std::endl;
    // std::cout << "max_num_face = " << max_num_face << std::endl;
    size_t num_iso_vert = 0;
    size_t num_iso_face = 0;
    //
    for (size_t i = 0; i < num_tet; i++) {
        if (has_isosurface[i]) {
            const auto& arrangement = cut_results[i];
            const auto& vertices = arrangement.vertices;
            const auto& faces = arrangement.faces;
            // const auto& func_ids = func_in_tet[i];
            //auto num_func = func_ids.size();  // num_func = 1
            // find vertices and faces on isosurface
            std::vector<bool> is_iso_vert(vertices.size(), false);
            std::vector<bool> is_iso_face(faces.size(), false);
            if (arrangement.unique_planes.empty()) { // all planes are unique
                for (size_t j = 0; j < faces.size(); ++j) {
                    if (faces[j].supporting_plane > 3) { // plane 0,1,2,3 are tet boundaries
                        is_iso_face[j] = true;
                        for (const auto& vid : faces[j].vertices) {
                            is_iso_vert[vid] = true;
                        }
                    }
                }
            } else {
                for (size_t j = 0; j < faces.size(); j++) {
                    auto pid = faces[j].supporting_plane;
                    auto uid = arrangement.unique_plane_indices[pid];
                    for (const auto& plane_id : arrangement.unique_planes[uid]) {
                        if (plane_id > 3) { // plane 0,1,2,3 are tet boundaries
                            is_iso_face[j] = true;
                            for (const auto& vid : faces[j].vertices) {
                                is_iso_vert[vid] = true;
                            }
                            break;
                        }
                    }
                }
            }
            // create iso-vertices
            std::vector<size_t> iso_vId_of_vert(vertices.size(), Arrangement<3>::None);
            for (size_t j = 0; j < vertices.size(); j++) {
                if (is_iso_vert[j]) {
                    // std::array<size_t, 3> implicit_pIds;
                    std::array<size_t, 3> bndry_pIds;
                    size_t num_bndry_planes = 0;
                    // size_t num_impl_planes = 0;
                    const auto& vertex = vertices[j];
                    // vertex.size() == 3
                    for (size_t k = 0; k < 3; k++) {
                        if (vertex[k] < 4) { // plane 0,1,2,3 are tet boundaries
                            bndry_pIds[num_bndry_planes] = vertex[k];
                            ++num_bndry_planes;
                        }
                    }
                    // assert: num_bndry_planes == 2 (on tet edge) or 3 (on tet vertex)
                    if (num_bndry_planes == 2) { // on tet edge
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
                        std::array<size_t, 2> key = {vId1, vId2};
                        auto iter_inserted = vert_on_tetEdge.try_emplace(key, num_iso_vert);
                        if (iter_inserted.second) {
                            auto& iso_vert = iso_verts[num_iso_vert];
                            iso_vert.tet_index = i;
                            iso_vert.tet_vert_index = j;
                            iso_vert.simplex_size = 2;
                            iso_vert.simplex_vert_indices[0] = vId1;
                            iso_vert.simplex_vert_indices[1] = vId2;
                            // iso_vert.func_indices[0] = implicit_pIds[0];
                            ++num_iso_vert;
                        }
                        iso_vId_of_vert[j] = iter_inserted.first->second;
                    } else { // on tet vertex
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
                    // face is on tet boundary if face.negative_cell is NONE
                    bool face_on_tet_boundary = (faces[j].negative_cell == Arrangement<3>::None);
                    if (face_on_tet_boundary) {
                        /*std::vector<size_t> sorted_face_verts = face_verts;
                        std::sort(sorted_face_verts.begin(), sorted_face_verts.end());
                        std::array<size_t, 3> key = {
                            sorted_face_verts[0], sorted_face_verts[1], sorted_face_verts[2]};*/
                        std::array<size_t, 3> key;
                        compute_iso_face_key(face_verts, key);
                        auto iter_inserted = face_on_tetFace.try_emplace(key, num_iso_face);
                        if (iter_inserted.second) {
                            iso_faces[num_iso_face].vert_indices = face_verts;
                            // iso_faces[num_iso_face].tet_face_indices.emplace_back(i, j);
                            ++num_iso_face;
                        } else { // iso_face inserted before
                            size_t iso_face_id = (iter_inserted.first)->second;
                            // iso_faces[iso_face_id].tet_face_indices.emplace_back(i, j);
                        }
                    } else { // face not on tet boundary
                        iso_faces[num_iso_face].vert_indices = face_verts;
                        // iso_faces[num_iso_face].tet_face_indices.emplace_back(i, j);
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

void compute_iso_vert_xyz(const std::vector<IsoVert>& iso_verts,
    const Eigen::MatrixXd& funcVals,
    const std::vector<std::array<double, 3>>& pts,
    std::vector<std::array<double, 3>>& iso_pts)
{
    iso_pts.resize(iso_verts.size());
    std::array<double, 2> b2;
    std::array<double, 3> f1s3;
    std::array<double, 3> f2s3;
    std::array<double, 3> b3;
    std::array<double, 4> f1s4;
    std::array<double, 4> f2s4;
    std::array<double, 4> f3s4;
    std::array<double, 4> b4;
    for (size_t i = 0; i < iso_verts.size(); i++) {
        const auto& iso_vert = iso_verts[i];
        switch (iso_vert.simplex_size) {
        case 2: // on tet edge
        {
            auto vId1 = iso_vert.simplex_vert_indices[0];
            auto vId2 = iso_vert.simplex_vert_indices[1];
            auto fId = iso_vert.func_indices[0];
            auto f1 = funcVals(fId, vId1);
            auto f2 = funcVals(fId, vId2);
            //
            compute_barycentric_coords(f1, f2, b2);
            iso_pts[i][0] = b2[0] * pts[vId1][0] + b2[1] * pts[vId2][0];
            iso_pts[i][1] = b2[0] * pts[vId1][1] + b2[1] * pts[vId2][1];
            iso_pts[i][2] = b2[0] * pts[vId1][2] + b2[1] * pts[vId2][2];
            break;
        }
        case 3: // on tet face
        {
            auto vId1 = iso_vert.simplex_vert_indices[0];
            auto vId2 = iso_vert.simplex_vert_indices[1];
            auto vId3 = iso_vert.simplex_vert_indices[2];
            auto fId1 = iso_vert.func_indices[0];
            auto fId2 = iso_vert.func_indices[1];
            //
            f1s3[0] = funcVals(fId1, vId1);
            f1s3[1] = funcVals(fId1, vId2);
            f1s3[2] = funcVals(fId1, vId3);
            //
            f2s3[0] = funcVals(fId2, vId1);
            f2s3[1] = funcVals(fId2, vId2);
            f2s3[2] = funcVals(fId2, vId3);
            //
            compute_barycentric_coords(f1s3, f2s3, b3);
            iso_pts[i][0] = b3[0] * pts[vId1][0] + b3[1] * pts[vId2][0] + b3[2] * pts[vId3][0];
            iso_pts[i][1] = b3[0] * pts[vId1][1] + b3[1] * pts[vId2][1] + b3[2] * pts[vId3][1];
            iso_pts[i][2] = b3[0] * pts[vId1][2] + b3[1] * pts[vId2][2] + b3[2] * pts[vId3][2];
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
            f1s4[0] = funcVals(fId1, vId1);
            f1s4[1] = funcVals(fId1, vId2);
            f1s4[2] = funcVals(fId1, vId3);
            f1s4[3] = funcVals(fId1, vId4);
            //
            f2s4[0] = funcVals(fId2, vId1);
            f2s4[1] = funcVals(fId2, vId2);
            f2s4[2] = funcVals(fId2, vId3);
            f2s4[3] = funcVals(fId2, vId4);
            //
            f3s4[0] = funcVals(fId3, vId1);
            f3s4[1] = funcVals(fId3, vId2);
            f3s4[2] = funcVals(fId3, vId3);
            f3s4[3] = funcVals(fId3, vId4);
            //
            compute_barycentric_coords(f1s4, f2s4, f3s4, b4);
            iso_pts[i][0] = b4[0] * pts[vId1][0] + b4[1] * pts[vId2][0] + b4[2] * pts[vId3][0] +
                            b4[3] * pts[vId4][0];
            iso_pts[i][1] = b4[0] * pts[vId1][1] + b4[1] * pts[vId2][1] + b4[2] * pts[vId3][1] +
                            b4[3] * pts[vId4][1];
            iso_pts[i][2] = b4[0] * pts[vId1][2] + b4[1] * pts[vId2][2] + b4[2] * pts[vId3][2] +
                            b4[3] * pts[vId4][2];
            break;
        }
        case 1: // on tet vertex
            iso_pts[i] = pts[iso_vert.simplex_vert_indices[0]];
            break;
        default: break;
        }
    }
}

void compute_iso_vert_xyz_marching_tet(const std::vector<IsoVert>& iso_verts,
    const std::vector<double>& funcVals,
    const std::vector<std::array<double, 3>>& pts,
    std::vector<std::array<double, 3>>& iso_pts)
{
    iso_pts.resize(iso_verts.size());

    for (size_t i = 0; i < iso_verts.size(); i++) {
        const auto& iso_vert = iso_verts[i];
        // assert: iso_vert.simplex_size == 2 (on tet edge) or 1 (on tet vertex)
        if (iso_vert.simplex_size == 2) { // on tet edge
            auto vId1 = iso_vert.simplex_vert_indices[0];
            auto vId2 = iso_vert.simplex_vert_indices[1];
            auto f1 = funcVals[vId1];
            auto f2 = funcVals[vId2];
            std::array<double, 2> b;
            compute_barycentric_coords(f1, f2, b);
            iso_pts[i][0] = b[0] * pts[vId1][0] + b[1] * pts[vId2][0];
            iso_pts[i][1] = b[0] * pts[vId1][1] + b[1] * pts[vId2][1];
            iso_pts[i][2] = b[0] * pts[vId1][2] + b[1] * pts[vId2][2];
        } else { // 1: on tet vertex
            iso_pts[i] = pts[iso_vert.simplex_vert_indices[0]];
        }
    }
}

//


// compute iso-edges and edge-face connectivity
// void compute_mesh_edges(std::vector<PolygonFace>& iso_faces, std::vector<Edge>& iso_edges)
//{
//    size_t max_num_edge = 0;
//    for (size_t i = 0; i < iso_faces.size(); i++) {
//        max_num_edge += iso_faces[i].vert_indices.size();
//    }
//    iso_edges.resize(max_num_edge);
//    size_t num_iso_edge = 0;
//    // map: (v1, v2) -> iso-edge index
//    absl::flat_hash_map<std::pair<size_t, size_t>, size_t> edge_id;
//    for (size_t i = 0; i < iso_faces.size(); i++) {
//        auto& face = iso_faces[i];
//        size_t num_edge = face.vert_indices.size();
//        face.edge_indices.resize(num_edge);
//        for (size_t j = 0; j < num_edge; j++) {
//            size_t v1 = face.vert_indices[j];
//            size_t v2 = (j + 1 == num_edge) ? face.vert_indices[0] : face.vert_indices[j + 1];
//            // swap if v1 > v2
//            size_t tmp = v1;
//            if (v1 > v2) {
//                v1 = v2;
//                v2 = tmp;
//            }
//            //
//            auto iter_inserted = edge_id.try_emplace(std::make_pair(v1, v2), num_iso_edge);
//            if (iter_inserted.second) { // new iso-edge
//                iso_edges[num_iso_edge].v1 = v1;
//                iso_edges[num_iso_edge].v2 = v2;
//                iso_edges[num_iso_edge].face_edge_indices.emplace_back(i, j);
//                face.edge_indices[j] = num_iso_edge;
//                ++num_iso_edge;
//            } else { // existing iso-edge
//                size_t eId = iter_inserted.first->second;
//                iso_edges[eId].face_edge_indices.emplace_back(i, j);
//                face.edge_indices[j] = eId;
//            }
//        }
//    }
//    //
//    iso_edges.resize(num_iso_edge);
//}

void compute_mesh_edges(const std::vector<PolygonFace>& mesh_faces,
    std::vector<std::vector<size_t>>& edges_of_face,
    std::vector<Edge>& mesh_edges)
{
    size_t max_num_edge = 0;
    for (const auto& iso_face : mesh_faces) {
        max_num_edge += iso_face.vert_indices.size();
    }
    mesh_edges.reserve(max_num_edge / 2);
    edges_of_face.reserve(mesh_faces.size());
    size_t num_iso_edge;
    // map: (v1, v2) -> iso-edge index
    absl::flat_hash_map<std::pair<size_t, size_t>, size_t> edge_id;
    for (size_t i = 0; i < mesh_faces.size(); i++) {
        auto& face = mesh_faces[i];
        size_t num_edge = face.vert_indices.size();
        //        face.edge_indices.resize(num_edge);
        edges_of_face.emplace_back(num_edge);
        auto& face_edges = edges_of_face.back();
        for (size_t j = 0; j < num_edge; j++) {
            size_t v1 = face.vert_indices[j];
            size_t v2 = (j + 1 == num_edge) ? face.vert_indices[0] : face.vert_indices[j + 1];
            // swap if v1 > v2
            size_t tmp = v1;
            if (v1 > v2) {
                v1 = v2;
                v2 = tmp;
            }
            //
            num_iso_edge = mesh_edges.size();
            auto iter_inserted = edge_id.try_emplace(std::make_pair(v1, v2), num_iso_edge);
            if (iter_inserted.second) { // new iso-edge
                mesh_edges.emplace_back();
                mesh_edges.back().v1 = v1;
                mesh_edges.back().v2 = v2;
                mesh_edges.back().face_edge_indices.emplace_back(i, j);
                //                face.edge_indices[j] = num_iso_edge;
                face_edges[j] = num_iso_edge;
            } else { // existing iso-edge
                size_t eId = iter_inserted.first->second;
                mesh_edges[eId].face_edge_indices.emplace_back(i, j);
                //                face.edge_indices[j] = eId;
                face_edges[j] = eId;
            }
        }
    }
}

// compute iso-edges and edge-face connectivity (reserve instead of resize)
// void compute_iso_edges_r(std::vector<PolygonFace>& iso_faces, std::vector<Edge>& iso_edges)
//{
//    size_t max_num_edge = 0;
//    for (size_t i = 0; i < iso_faces.size(); i++) {
//        max_num_edge += iso_faces[i].vert_indices.size();
//    }
//    iso_edges.reserve(max_num_edge/2);
//    size_t num_iso_edge = 0;
//    // map: (v1, v2) -> iso-edge index
//    absl::flat_hash_map<std::pair<size_t, size_t>, size_t> edge_id;
//    for (size_t i = 0; i < iso_faces.size(); i++) {
//        auto& face = iso_faces[i];
//        size_t num_edge = face.vert_indices.size();
//        face.edge_indices.resize(num_edge);
//        for (size_t j = 0; j < num_edge; j++) {
//            size_t v1 = face.vert_indices[j];
//            size_t v2 = (j + 1 == num_edge) ? face.vert_indices[0] : face.vert_indices[j + 1];
//            // swap if v1 > v2
//            size_t tmp = v1;
//            if (v1 > v2) {
//                v1 = v2;
//                v2 = tmp;
//            }
//            //
//            num_iso_edge = iso_edges.size();
//            auto iter_inserted = edge_id.try_emplace(std::make_pair(v1, v2), num_iso_edge);
//            if (iter_inserted.second) { // new iso-edge
//                iso_edges.emplace_back();
//                iso_edges.back().v1 = v1;
//                iso_edges.back().v2 = v2;
//                iso_edges.back().face_edge_indices.emplace_back(i,j);
//                face.edge_indices[j] = num_iso_edge;
//            } else { // existing iso-edge
//                size_t eId = iter_inserted.first->second;
//                iso_edges[eId].face_edge_indices.emplace_back(i, j);
//                face.edge_indices[j] = eId;
//            }
//        }
//    }
//}

// group iso-faces into patches
// void compute_patches(const std::vector<PolygonFace>& iso_faces,
//    const std::vector<Edge>& iso_edges,
//    std::vector<std::vector<size_t>>& patches)
//{
//    std::vector<bool> visisted_face(iso_faces.size(), false);
//    for (size_t i = 0; i < iso_faces.size(); i++) {
//        if (!visisted_face[i]) {
//            // new patch
//            patches.emplace_back();
//            auto& patch = patches.back();
//            std::queue<size_t> Q;
//            Q.push(i);
//            patch.emplace_back(i);
//            visisted_face[i] = true;
//            while (!Q.empty()) {
//                auto fId = Q.front();
//                Q.pop();
//                for (size_t eId : iso_faces[fId].edge_indices) {
//                    if (iso_edges[eId].face_edge_indices.size() == 2) { // manifold edge
//                        size_t other_fId = (iso_edges[eId].face_edge_indices[0].first == fId)
//                                               ? iso_edges[eId].face_edge_indices[1].first
//                                               : iso_edges[eId].face_edge_indices[0].first;
//                        if (!visisted_face[other_fId]) {
//                            Q.push(other_fId);
//                            patch.emplace_back(other_fId);
//                            visisted_face[other_fId] = true;
//                        }
//                    }
//                }
//            }
//        }
//    }
//}

void compute_patches(const std::vector<std::vector<size_t>>& edges_of_face,
    const std::vector<Edge>& mesh_edges,
    std::vector<std::vector<size_t>>& patches)
{
    std::vector<bool> visisted_face(edges_of_face.size(), false);
    for (size_t i = 0; i < edges_of_face.size(); i++) {
        if (!visisted_face[i]) {
            // new patch
            patches.emplace_back();
            auto& patch = patches.back();
            std::queue<size_t> Q;
            Q.push(i);
            patch.emplace_back(i);
            visisted_face[i] = true;
            while (!Q.empty()) {
                auto fId = Q.front();
                Q.pop();
                for (size_t eId : edges_of_face[fId]) {
                    if (mesh_edges[eId].face_edge_indices.size() == 2) { // manifold edge
                        size_t other_fId = (mesh_edges[eId].face_edge_indices[0].first == fId)
                                               ? mesh_edges[eId].face_edge_indices[1].first
                                               : mesh_edges[eId].face_edge_indices[0].first;
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
void compute_chains(const std::vector<Edge>& mesh_edges,
    const std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
    std::vector<std::vector<size_t>>& chains)
{
    std::vector<bool> visited_edge(mesh_edges.size(), false);
    for (size_t i = 0; i < mesh_edges.size(); i++) {
        if (!visited_edge[i] && mesh_edges[i].face_edge_indices.size() > 2) {
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
                size_t v = mesh_edges[eId].v1;
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
                v = mesh_edges[eId].v2;
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
    const std::vector<PolygonFace>& iso_faces,
    const Edge& iso_edge,
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

void compute_shells_and_components(size_t num_patch,
    const std::vector<std::vector<std::pair<size_t, int>>>& half_patch_list,
    std::vector<std::vector<size_t>>& shells,
    std::vector<size_t>& shell_of_half_patch,
    std::vector<std::vector<size_t>>& components,
    std::vector<size_t>& component_of_patch)
{
    // (patch i, 1) <--> 2i,  (patch i, -1) <--> 2i+1
    // compute half-patch adjacency list
    std::vector<std::vector<size_t>> half_patch_adj_list(2 * num_patch);
    for (const auto& half_patches : half_patch_list) {
        for (size_t i = 0; i + 1 < half_patches.size(); i++) {
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
    // each component is a shell
    std::vector<bool> visited_half_patch(2 * num_patch, false);
    shells.clear();
    shell_of_half_patch.resize(2 * num_patch);
    for (size_t i = 0; i < 2 * num_patch; i++) {
        if (!visited_half_patch[i]) {
            // create new component
            size_t shell_Id = shells.size();
            shells.emplace_back();
            auto& shell = shells.back();
            std::queue<size_t> Q;
            Q.push(i);
            shell.push_back(i);
            shell_of_half_patch[i] = shell_Id;
            visited_half_patch[i] = true;
            while (!Q.empty()) {
                auto half_patch = Q.front();
                Q.pop();
                for (auto hp : half_patch_adj_list[half_patch]) {
                    if (!visited_half_patch[hp]) {
                        shell.push_back(hp);
                        shell_of_half_patch[hp] = shell_Id;
                        Q.push(hp);
                        visited_half_patch[hp] = true;
                    }
                }
            }
        }
    }
    // find connected component of patch-adjacency graph
    // each component is a iso-surface component
    std::vector<bool> visited_patch(num_patch, false);
    components.clear();
    component_of_patch.resize(num_patch);
    for (size_t i = 0; i < num_patch; i++) {
        if (!visited_patch[i]) {
            // create new component
            size_t component_Id = components.size();
            components.emplace_back();
            auto& component = components.back();
            std::queue<size_t> Q;
            Q.push(i);
            component.push_back(i);
            component_of_patch[i] = component_Id;
            visited_patch[i] = true;
            while (!Q.empty()) {
                auto patch = Q.front();
                Q.pop();
                // 1-side of patch
                for (auto hp : half_patch_adj_list[2 * patch]) {
                    if (!visited_patch[hp/2]) {
                        auto p = hp/2;
                        component.push_back(p);
                        component_of_patch[p] = component_Id;
                        Q.push(p);
                        visited_patch[p] = true;
                    }
                }
                // -1-side of patch
                for (auto hp : half_patch_adj_list[2 * patch + 1]) {
                    if (!visited_patch[hp/2]) {
                        auto p = hp/2;
                        component.push_back(p);
                        component_of_patch[p] = component_Id;
                        Q.push(p);
                        visited_patch[p] = true;
                    }
                }
            }
        }
    }
    // get shells as list of patch indices
//    for (auto& shell : shells) {
//        for (auto& pId : shell) {
//            pId /= 2; // get patch index of half-patch
//        }
//    }
}


void compute_edge_intersection_order(
    const Arrangement<3>& tet_cut_result, size_t v, size_t u, std::vector<size_t>& vert_indices)
{
    const auto& vertices = tet_cut_result.vertices;
    const auto& faces = tet_cut_result.faces;
    //
    std::array<bool,4> edge_flag {true, true, true, true};
    edge_flag[v] = false;
    edge_flag[u] = false;

    // find vertices on edge v->u, and index of v and u in tet_cut_result.vertices
    size_t v_id, u_id;
    std::vector<bool> is_on_edge_vu(vertices.size(), false);
    size_t num_vert_in_edge_vu = 0;
    size_t flag_count;
    size_t other_plane;
    for (size_t i = 0; i < vertices.size(); ++i) {
        flag_count = 0;
        other_plane = Arrangement<3>::None;
        const auto& vert = vertices[i];
        for (int j = 0; j < 3; ++j) {
            if (vert[j] < 4) {  // 0,1,2,3 are tet boundary planes
                if (edge_flag[vert[j]]) {
                    ++flag_count;
                } else {
                   other_plane = vert[j];
                }
            }
        }
        if (flag_count == 2) {
            is_on_edge_vu[i] = true;
            if (other_plane == u) {  // current vertex is v
                v_id = i;
            } else if (other_plane == v) {  // current vertex is u
                u_id = i;
            } else {  // current vertex is in interior of edge v->u
                ++num_vert_in_edge_vu;
            }
        }
    }
    if (num_vert_in_edge_vu == 0) { // no intersection points in edge v->u
        vert_indices.push_back(v_id);
        vert_indices.push_back(u_id);
        return;
    }

    // find all faces on triangle v->u->w, w is a third tet vertex
    size_t pId;  // plane index of v->u->w, pick pId != v and pId != u
    for (size_t i = 0; i < 4; ++i) {
        if (edge_flag[i]) {
            pId = i;
            break;
        }
    }
//    std::vector<bool> is_on_tri_vuw(faces.size(), false);
    std::vector<size_t> faces_on_tri_vuw;
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& face = faces[i];
        if (face.negative_cell == Arrangement<3>::None) { // face is on tet boundary
            if (face.supporting_plane == pId) {
                faces_on_tri_vuw.push_back(i);
            }
        }
    }

    // build edge-face connectivity on triangle v->u->w
    // map: edge (v1, v2) --> incident faces (f1, f2), f2 is None if (v1,v2) is on triangle boundary
    absl::flat_hash_map<std::pair<size_t, size_t>, std::pair<size_t,size_t>> faces_of_edge;
    for (auto fId : faces_on_tri_vuw) {
        const auto& face = faces[fId];
        size_t num_vert = face.vertices.size();
        for (size_t i = 0; i < num_vert; ++i) {
            size_t i_next = (i + 1) % num_vert;
            size_t vi = face.vertices[i];
            size_t vi_next = face.vertices[i_next];
            // add fId to edge (vi, vi_next)
            auto iter_inserted = faces_of_edge.try_emplace(
                std::make_pair(vi, vi_next),
                std::make_pair(fId, Arrangement<3>::None));
            if (!iter_inserted.second) {  // inserted before
                iter_inserted.first->second.second = fId;
            }
            // add fId to edge (vi_next, vi)
            iter_inserted = faces_of_edge.try_emplace(
                std::make_pair(vi_next, vi),
                std::make_pair(fId, Arrangement<3>::None));
            if (!iter_inserted.second) {  // inserted before
                iter_inserted.first->second.second = fId;
            }
        }
    }

    // find the face on triangle v->u->w:
    // 1. has vertex v
    // 2. has an edge on v->u
    size_t f_start;
    for (auto fId : faces_on_tri_vuw) {
        const auto& face = faces[fId];
        bool find_v = false;
        size_t count = 0;
        for (auto vi : face.vertices) {
            if (vi == v_id) {
                find_v = true;
            }
            if (is_on_edge_vu[vi]) {
                ++count;
            }
        }
        if (find_v && count == 2) {
            f_start = fId;
            break;
        }
    }

    // trace edge v->u
    vert_indices.reserve(num_vert_in_edge_vu + 2);
    vert_indices.push_back(v_id);
    std::vector<bool> visited_face(faces.size(), false);
    size_t v_curr = v_id;
    size_t f_curr = f_start;
    std::pair<size_t, size_t> edge_prev;
    std::pair<size_t, size_t> edge_next;
    std::pair<size_t, size_t> edge_on_vu;
    while (v_curr != u_id) {
        // clear edge_on_vu
        edge_on_vu.first = Arrangement<3>::None;
        edge_on_vu.second = Arrangement<3>::None;
        //
        const auto& face = faces[f_curr];
        size_t num_vert = face.vertices.size();
        // visit all edges of face, find edge_prev, edge_next and edge_on_vu
        for (size_t i = 0; i < num_vert; ++i) {
            size_t i_next = (i+1) % num_vert;
            size_t vi = face.vertices[i];
            size_t vi_next = face.vertices[i_next];
            if (is_on_edge_vu[vi] && !is_on_edge_vu[vi_next]) {
                auto& two_faces = faces_of_edge[std::make_pair(vi, vi_next)];
                auto other_face = (two_faces.first == f_curr) ? two_faces.second : two_faces.first;
                if (vi == v_id || (other_face != Arrangement<3>::None && visited_face[other_face])) {
                    edge_prev.first = vi;
                    edge_prev.second = vi_next;
                } else {
                    edge_next.first = vi;
                    edge_next.second = vi_next;
                }
            } else if (is_on_edge_vu[vi_next] && !is_on_edge_vu[vi]) {
                auto& two_faces = faces_of_edge[std::make_pair(vi, vi_next)];
                auto other_face = (two_faces.first == f_curr) ? two_faces.second : two_faces.first;
                if (vi_next == v_id || (other_face != Arrangement<3>::None && visited_face[other_face])) {
                    edge_prev.first = vi;
                    edge_prev.second = vi_next;
                } else {
                    edge_next.first = vi;
                    edge_next.second = vi_next;
                }
            } else if (is_on_edge_vu[vi] && is_on_edge_vu[vi_next]) {
                edge_on_vu.first = vi;
                edge_on_vu.second = vi_next;
            }
        }
        //
        if (edge_on_vu.first == Arrangement<3>::None) {
            // no edge of the face is on v->u
            // keep v_curr, update f_curr
            visited_face[f_curr] = true;
            auto& two_faces = faces_of_edge[edge_next];
            f_curr = (two_faces.first == f_curr) ? two_faces.second : two_faces.first;
        } else {
            // there is an edge of the face on v->u
            // update v_curr to be the next vert on edge v->u
            v_curr = (edge_on_vu.first == v_curr) ? edge_on_vu.second : edge_on_vu.first;
            vert_indices.push_back(v_curr);
            // update f_curr
            visited_face[f_curr] = true;
            auto& two_faces = faces_of_edge[edge_next];
            f_curr = (two_faces.first == f_curr) ? two_faces.second : two_faces.first;
        }
    }
}


void compute_passing_face(
    const Arrangement<3>& tet_cut_result, size_t v, size_t u, std::pair<size_t, int>& face_orient){
    // find a face incident to edge v -> u
    const auto& faces = tet_cut_result.faces;
    size_t incident_face_id;
    bool found_incident_face = false;
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& face = faces[i];
        size_t num_vert = face.vertices.size();
        for (size_t j = 0; j < num_vert; ++j) {
            size_t vId1 = face.vertices[j];
            size_t vId2 = face.vertices[(j+1) % num_vert];
            if ((vId1==v && vId2==u) || (vId1==u && vId2 == v)) {
                incident_face_id = i;
                found_incident_face = true;
                break;
            }
        }
        if (found_incident_face) {
            break;
        }
    }
    // assert: found_incident_face == true
    size_t cell_id = faces[incident_face_id].positive_cell;
    const auto& cell = tet_cut_result.cells[cell_id];

    // find the face
    // 1. bounding the cell
    // 2. passing vertex v
    for (auto fId : cell.faces) {
        const auto& face = faces[fId];
        bool found_v = false;
        bool found_u = false;
        for (auto vId : face.vertices) {
            if (vId == v) {
                found_v = true;
            } else if (vId == u) {
                found_u = true;
            }
        }
        if (found_v && !found_u) {
            // the current face passes v but not u
            face_orient.first = fId;
            face_orient.second = (face.positive_cell == cell_id) ? 1 : -1;
            return;
        }
    }
}

void compute_passing_face_pair(const Arrangement<3>& tet_cut_result,
    size_t v1,
    size_t v2,
    std::pair<size_t, int>& face_orient1,
    std::pair<size_t, int>& face_orient2){
    // find a face incident to edge v1 -> v2
    const auto& faces = tet_cut_result.faces;
    size_t incident_face_id;
    bool found_incident_face = false;
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& face = faces[i];
        size_t num_vert = face.vertices.size();
        for (size_t j = 0; j < num_vert; ++j) {
            size_t vId1 = face.vertices[j];
            size_t vId2 = face.vertices[(j+1) % num_vert];
            if ((vId1==v1 && vId2==v2) || (vId1==v2 && vId2 == v1)) {
                incident_face_id = i;
                found_incident_face = true;
                break;
            }
        }
        if (found_incident_face) {
            break;
        }
    }
    // assert: found_incident_face == true
    size_t cell_id = faces[incident_face_id].positive_cell;
    const auto& cell = tet_cut_result.cells[cell_id];

    // find the two faces
    // 1. bounding the cell
    // 2. passing vertex v1 or v2
    for (auto fId : cell.faces) {
        const auto& face = faces[fId];
        bool found_v1 = false;
        bool found_v2 = false;
        for (auto vId : face.vertices) {
            if (vId == v1) {
                found_v1 = true;
            } else if (vId == v2) {
                found_v2 = true;
            }
        }
        if (found_v1 && !found_v2) {
            // the current face passes v1 but not v2
            face_orient1.first = fId;
            face_orient1.second = (face.positive_cell == cell_id) ? 1 : -1;
        } else if (!found_v1 && found_v2) {
            // the current face passes v2 but not v1
            face_orient2.first = fId;
            face_orient2.second = (face.positive_cell == cell_id) ? 1 : -1;
        }
    }
}


void compute_arrangement_cells(size_t num_shell,
    const std::vector<std::pair<size_t, size_t>>& shell_links,
    std::vector<std::vector<size_t>>& arrangement_cells){
    // build shell adjacency list
    size_t sink_shell = num_shell;
    absl::flat_hash_map<size_t, std::vector<size_t>> adjacent_shells;
    for (const auto& link : shell_links) {
        if (link.first == Arrangement<3>::None) {
            adjacent_shells[sink_shell].push_back(link.second);
            adjacent_shells[link.second].push_back(sink_shell);
        } else if (link.second == Arrangement<3>::None) {
            adjacent_shells[sink_shell].push_back(link.first);
            adjacent_shells[link.first].push_back(sink_shell);
        } else {
            adjacent_shells[link.first].push_back(link.second);
            adjacent_shells[link.second].push_back(link.first);
        }
    }

    // find connected components of shell adjacency graph
    // each component is an arrangement cells
    std::vector<bool> visited_shell(num_shell+1, false);
//    arrangement_cells.clear();
    for (size_t i = 0; i < num_shell+1 ; ++i) {
        if (!visited_shell[i]) {
            // create new component
            arrangement_cells.emplace_back();
            auto& arr_cell = arrangement_cells.back();
            std::queue<size_t> Q;
            Q.push(i);
            arr_cell.push_back(i);
            visited_shell[i] = true;
            while (!Q.empty()) {
                auto shell_id = Q.front();
                Q.pop();
                for (auto s : adjacent_shells[shell_id]) {
                    if (!visited_shell[s]) {
                        arr_cell.push_back(s);
                        Q.push(s);
                        visited_shell[s] = true;
                    }
                }
            }
        }
    }

    // remove sink shell from arrangement cells
    std::vector<size_t> sink_free_shell_list;
    for (auto& arr_cell : arrangement_cells) {
        sink_free_shell_list.clear();
        for (auto s : arr_cell) {
            if (s < num_shell) {
                sink_free_shell_list.push_back(s);
            }
        }
        arr_cell = sink_free_shell_list;
    }
}
bool parse_config_file_MI(const std::string& filename,
    std::string& tet_mesh_file,
    std::string& material_file,
    std::string& output_dir,
    bool& use_lookup,
    bool& use_3func_lookup)
{
    using json = nlohmann::json;
    std::ifstream fin(filename.c_str());
    if (!fin) {
        std::cout << "configure file not exist!" << std::endl;
        return false;
    }
    json data;
    fin >> data;
    fin.close();
    //
    tet_mesh_file = data["tetMeshFile"];
    material_file = data["materialFile"];
    output_dir = data["outputDir"];
    use_lookup = data["useLookup"];
    use_3func_lookup = data["use3funcLookup"];
    return true;
}

bool load_seeds(const std::string& filename, std::vector<std::array<double, 3>>& seeds)
{
    using json = nlohmann::json;
    std::ifstream fin(filename.c_str());
    if (!fin) {
        std::cout << "seed file not exist!" << std::endl;
        return false;
    }
    json data;
    fin >> data;
    fin.close();
    //
    seeds.resize(data.size());
    for (size_t i = 0; i < seeds.size(); i++) {
        // seed: {x,y,z}
        for (int j = 0; j < 3; ++j) {
            seeds[i][j] = data[i][j].get<double>();
        }
    }
    return true;
};


void extract_MI_mesh_pure(size_t num_2_func,
    size_t num_3_func,
    size_t num_more_func,
    const std::vector<MaterialInterface<3>>& cut_results,
    const std::vector<size_t>& cut_result_index,
    const std::vector<size_t>& material_in_tet,
    const std::vector<size_t>& start_index_of_tet,
    const std::vector<std::array<size_t, 4>>& tets,
    std::vector<MI_Vert>& MI_verts,
    std::vector<PolygonFace>& MI_faces)
{
    size_t n_tets = tets.size();
    // estimate number of MI-verts and MI-faces
    // C(2,2) = 1, C(3,2) = 3, C(4,2) = 6
    size_t max_num_face = num_2_func + 3 * num_3_func + 6 * num_more_func;
    // for triangle mesh, 3|F| is about 2|F|, |V|-|E|+|F| is about 0, so |V| is about 0.5|F|.
    // Since most non-empty tets has 2 functions (i.e. 1 material interface),
    // most polygon's in the mesh are triangle or quad.
    // so number of triangles is roughly between |P| and 2|P| (|P| = number of polygons)
    // so, 0.5 * 2|P| = |P| is an estimate of |V|'s upper-bound
    size_t max_num_vert = max_num_face;
    MI_verts.reserve(max_num_vert);
    MI_faces.reserve(max_num_face);
    // hash table for vertices on the boundary of tetrahedron
    absl::flat_hash_map<size_t, size_t> vert_on_tetVert;
    // key is (edge v1, edge v2, material 1, material 2)
    absl::flat_hash_map<std::array<size_t, 4>, size_t> vert_on_tetEdge;
    vert_on_tetEdge.reserve(num_2_func + 2 * num_3_func + 4 * num_more_func);
    // key is (tri v1, tri v2, tri v3, material 1, material 2, material 3)
    absl::flat_hash_map<std::array<size_t, 6>, size_t> vert_on_tetFace;
    // 2func: 0 face vert, 3func: 2 face vert, 4func: 4 face vert
    // num face vert estimate = (2 * n_3func + 4 * n_moreFunc) /2
    vert_on_tetFace.reserve(num_3_func + 2 * num_more_func);
    // map: face on tet boundary -> material label
    absl::flat_hash_map<std::array<long long,3>, size_t> material_of_bndry_face;
    //
    std::vector<bool> is_MI_vert;
    is_MI_vert.reserve(8);
    std::vector<bool> is_MI_face;
    is_MI_face.reserve(9);
    std::vector<long long> MI_vId_of_vert;
    MI_vId_of_vert.reserve(8);
    std::vector<size_t> face_verts;
    face_verts.reserve(4);
    std::vector<long long> bndry_face_verts;
    bndry_face_verts.reserve(4);
    std::array<long long,3> key3;
    std::array<size_t, 4> key4;
    std::array<size_t, 6> key6;
    std::array<bool, 4> used_mId;
    std::array<size_t, 2> vIds2;
    std::array<size_t, 3> vIds3;
    std::array<size_t, 4> func_mIds;
    std::array<size_t, 3> bndry_mIds;
    //
    for (size_t i = 0; i < n_tets; i++) {
        if (cut_result_index[i] != MaterialInterface<3>::None) {
            const auto& mInterface = cut_results[cut_result_index[i]];
            const auto& vertices = mInterface.vertices;
            const auto& faces = mInterface.faces;
            auto start_index = start_index_of_tet[i];
            auto num_func = start_index_of_tet[i + 1] - start_index;
            // find vertices and faces on material interface
            is_MI_vert.clear();
            for (int j = 0; j < vertices.size(); ++j) {
                is_MI_vert.push_back(false);
            }
            is_MI_face.clear();
            for (const auto & face : faces) {
                // material label 0,1,2,3 represents tet boundary
                if (face.positive_material_label > 3) {
                    is_MI_face.push_back(true);
                    for (auto vId : face.vertices) {
                        is_MI_vert[vId] = true;
                    }
                } else {
                    is_MI_face.push_back(false);
                }
            }
            // map: local vert index --> MI-vert index
            MI_vId_of_vert.clear();
            // create MI-vertices
            for (size_t j = 0; j < vertices.size(); j++) {
                size_t num_bndry_materials = 0;
                size_t num_func_materials = 0;
                const auto& vertex = vertices[j];
                // vertex.size() == 4
                for (size_t k = 0; k < 4; k++) {
                    if (vertex[k] > 3) { // material 0,1,2,3 are tet boundaries
                        func_mIds[num_func_materials] =
                            material_in_tet[vertex[k] - 4 + start_index];
                        ++num_func_materials;
                    } else {
                        bndry_mIds[num_bndry_materials] = vertex[k];
                        ++num_bndry_materials;
                    }
                }
                if (!is_MI_vert[j]) {
                    // vert not on any interior face, therefore must be a tet vertex
                    used_mId[0] = false;
                    used_mId[1] = false;
                    used_mId[2] = false;
                    used_mId[3] = false;
                    for (const auto& mId : bndry_mIds) {
                        used_mId[mId] = true;
                    }
                    size_t vId;
                    for (size_t k = 0; k < 4; k++) {
                        if (!used_mId[k]) {
                            vId = k;
                            break;
                        }
                    }
                    // tet vert i is mapped to (-i-1)
                    MI_vId_of_vert.push_back(-tets[i][vId] - 1);
                }
                else { // vert on an interior face
                    switch (num_bndry_materials) {
                    case 2: // on tet edge
                    {
                        used_mId[0] = false;
                        used_mId[1] = false;
                        used_mId[2] = false;
                        used_mId[3] = false;
                        used_mId[bndry_mIds[0]] = true;
                        used_mId[bndry_mIds[1]] = true;
                        size_t num_vIds = 0;
                        for (size_t k = 0; k < 4; k++) {
                            if (!used_mId[k]) {
                                vIds2[num_vIds] = tets[i][k];
                                ++num_vIds;
                            }
                        }
                        // sort {vId1, vId2}
                        size_t vId1 = vIds2[0];
                        size_t vId2 = vIds2[1];
                        if (vId1 > vId2) {
                            size_t tmp = vId1;
                            vId1 = vId2;
                            vId2 = tmp;
                        }
                        // sort {mId1, mId2}
                        size_t mId1 = func_mIds[0];
                        size_t mId2 = func_mIds[1];
                        if (mId1 > mId2) {
                            size_t tmp = mId1;
                            mId1 = mId2;
                            mId2 = tmp;
                        }
                        key4[0] = vId1;
                        key4[1] = vId2;
                        key4[2] = mId1;
                        key4[3] = mId2;
                        auto iter_inserted = vert_on_tetEdge.try_emplace(key4, MI_verts.size());
                        if (iter_inserted.second) {
                            MI_verts.emplace_back();
                            auto& MI_vert = MI_verts.back();
                            MI_vert.tet_index = i;
                            MI_vert.tet_vert_index = j;
                            MI_vert.simplex_size = 2;
                            MI_vert.simplex_vert_indices[0] = vId1;
                            MI_vert.simplex_vert_indices[1] = vId2;
                            MI_vert.material_indices[0] = mId1;
                            MI_vert.material_indices[1] = mId2;
                        }
                        MI_vId_of_vert.push_back(iter_inserted.first->second);
                        break;
                    }
                    case 1: // on tet face
                    {
                        size_t mId = bndry_mIds[0];
                        size_t num_vIds = 0;
                        for (size_t k = 0; k < 4; k++) {
                            if (k != mId) {
                                vIds3[num_vIds] = tets[i][k];
                                ++num_vIds;
                            }
                        }
                        std::sort(vIds3.begin(), vIds3.end());
                        std::sort(func_mIds.begin(), func_mIds.begin() + 3);
                        key6[0] = vIds3[0];
                        key6[1] = vIds3[1];
                        key6[2] = vIds3[2];
                        key6[3] = func_mIds[0];
                        key6[4] = func_mIds[1];
                        key6[5] = func_mIds[2];
                        auto iter_inserted = vert_on_tetFace.try_emplace(key6, MI_verts.size());
                        if (iter_inserted.second) {
                            MI_verts.emplace_back();
                            auto& MI_vert = MI_verts.back();
                            MI_vert.tet_index = i;
                            MI_vert.tet_vert_index = j;
                            MI_vert.simplex_size = 3;
                            MI_vert.simplex_vert_indices[0] = vIds3[0];
                            MI_vert.simplex_vert_indices[1] = vIds3[1];
                            MI_vert.simplex_vert_indices[2] = vIds3[2];
                            MI_vert.material_indices[0] = func_mIds[0];
                            MI_vert.material_indices[1] = func_mIds[1];
                            MI_vert.material_indices[2] = func_mIds[2];
                        }
                        MI_vId_of_vert.push_back(iter_inserted.first->second);
                        break;
                    }
                    case 0: // in tet cell
                    {
                        MI_vId_of_vert.push_back(MI_verts.size());
                        MI_verts.emplace_back();
                        auto& MI_vert = MI_verts.back();
                        MI_vert.tet_index = i;
                        MI_vert.tet_vert_index = j;
                        MI_vert.simplex_size = 4;
                        MI_vert.simplex_vert_indices = tets[i];
                        MI_vert.material_indices = func_mIds;
                        break;
                    }
                    case 3: // on tet vertex
                    {
                        used_mId[0] = false;
                        used_mId[1] = false;
                        used_mId[2] = false;
                        used_mId[3] = false;
                        for (const auto& mId : bndry_mIds) {
                            used_mId[mId] = true;
                        }
                        size_t vId;
                        for (size_t k = 0; k < 4; k++) {
                            if (!used_mId[k]) {
                                vId = k;
                                break;
                            }
                        }
                        auto key = tets[i][vId];
                        auto iter_inserted = vert_on_tetVert.try_emplace(key, MI_verts.size());
                        if (iter_inserted.second) {
                            MI_verts.emplace_back();
                            auto& MI_vert = MI_verts.back();
                            MI_vert.tet_index = i;
                            MI_vert.tet_vert_index = j;
                            MI_vert.simplex_size = 1;
                            MI_vert.simplex_vert_indices[0] = tets[i][vId];
                        }
                        MI_vId_of_vert.push_back(iter_inserted.first->second);
                        break;
                    }
                    default: break;
                    }
                }
            }
            // create MI-faces
            for (size_t j = 0; j < faces.size(); j++) {
                const auto& face = faces[j];
                if (is_MI_face[j]) {  // face j is in tet interior
                    face_verts.clear();
                    for (unsigned long vId : face.vertices) {
                        face_verts.push_back(MI_vId_of_vert[vId]);
                    }
                    MI_faces.emplace_back();
                    MI_faces.back().vert_indices = face_verts;
                    MI_faces.back().tet_face_indices.emplace_back(i, j);
                } else {  // face j is on tet boundary
                    bndry_face_verts.clear();
                    for (unsigned long vId : face.vertices) {
                        if (MI_vId_of_vert[vId] >= 0 && MI_verts[MI_vId_of_vert[vId]].simplex_size == 1) {
                            // intersection on tet vertex
                            bndry_face_verts.push_back(-MI_verts[MI_vId_of_vert[vId]].simplex_vert_indices[0] - 1);
                        } else {
                            bndry_face_verts.push_back(MI_vId_of_vert[vId]);
                        }
                    }
                    compute_iso_face_key(bndry_face_verts, key3);
                    auto m = material_in_tet[face.negative_material_label - 4 + start_index];
                    auto iter_inserted = material_of_bndry_face.try_emplace(key3, m);
                    if (! iter_inserted.second) { // inserted before
                        if (iter_inserted.first->second == m) {
                            // material labels on both sides of the face are the same
                            // the face is not on material interface
                            material_of_bndry_face.erase(iter_inserted.first);
                        } else {
                            // different material labels on different sides
                            // the face is on material interface
                            face_verts.clear();
                            for (size_t k = 0; k < face.vertices.size(); ++k) {
                                if (bndry_face_verts[k] < 0) {
                                        // try to create new vert on tet vertex
                                        size_t key = -bndry_face_verts[k] - 1;
                                        auto iter_inserted =
                                            vert_on_tetVert.try_emplace(key, MI_verts.size());
                                        if (iter_inserted.second) {
                                            MI_verts.emplace_back();
                                            auto& MI_vert = MI_verts.back();
                                            MI_vert.tet_index = i;
                                            MI_vert.tet_vert_index = j;
                                            MI_vert.simplex_size = 1;
                                            MI_vert.simplex_vert_indices[0] = key;
                                        }
                                        //                                    MI_vId_of_vert.push_back(iter_inserted.first->second);
                                    face_verts.push_back(iter_inserted.first->second);
                                } else {
                                    face_verts.push_back(MI_vId_of_vert[face.vertices[k]]);
                                }
                            }
                            // insert the face to MI_faces
                            MI_faces.emplace_back();
                            MI_faces.back().vert_indices = face_verts;
                            MI_faces.back().tet_face_indices.emplace_back(i, j);
                        }
                    }
                }
            }
        }
    }
}

bool save_result_MI(const std::string& filename,
    const std::vector<std::array<double, 3>>& MI_pts,
    const std::vector<PolygonFace>& MI_faces,
    const std::vector<std::vector<size_t>>& patches,
    const std::vector<Edge>& edges,
    const std::vector<std::vector<size_t>>& chains,
    const std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
    const std::vector<std::vector<std::pair<size_t, int>>>& half_patch_list,
    const std::vector<std::vector<size_t>>& shells,
    const std::vector<std::vector<size_t>>& components,
    const std::vector<std::vector<size_t>>& material_cells)
{
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str());
    //
    json jPts;
    for (const auto& pt : MI_pts) {
        jPts.push_back(json(pt));
    }
    //
    json jFaces;
    for (const auto& face : MI_faces) {
        jFaces.push_back(json(face.vert_indices));
    }
    //
    json jPatches;
    for (const auto& patch : patches) {
        jPatches.push_back(json(patch));
    }
    //
    json jEdges;
    for (const auto& edge : edges) {
        jEdges.push_back({edge.v1, edge.v2});
    }
    //
    json jChains;
    for (const auto& chain : chains) {
        jChains.push_back(json(chain));
    }
    //
    json jCorners;
    for (size_t i = 0; i < non_manifold_edges_of_vert.size(); i++) {
        if (non_manifold_edges_of_vert[i].size() > 2 || non_manifold_edges_of_vert[i].size() == 1) {
            jCorners.push_back(i);
        }
    }
    //
    json jHalfPatchsList;
    for (const auto& i : half_patch_list) {
        json jHalfPatchs;
        for (const auto& j : i) {
            jHalfPatchs.push_back(json(j));
        }
        jHalfPatchsList.push_back(jHalfPatchs);
    }
    //
    json jShells;
    for (const auto & shell : shells) {
        jShells.push_back(json(shell));
    }
    //
    json jComponents;
    for (const auto & component : components) {
        jComponents.push_back(json(component));
    }
    //
    json jMatCells;
    for (const auto& material_cell : material_cells) {
        jMatCells.push_back(json(material_cell));
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
    jOut.push_back(jShells);
    jOut.push_back(jComponents);
    jOut.push_back(jMatCells);
    fout << jOut << std::endl;
    fout.close();
    return true;
}

void compute_MI_vert_xyz(const std::vector<MI_Vert>& MI_verts,
    const Eigen::MatrixXd& funcVals,
    const std::vector<std::array<double, 3>>& pts,
    std::vector<std::array<double, 3>>& MI_pts)
{
    MI_pts.resize(MI_verts.size());
    std::array<double, 2> b2;
    std::array<double, 3> f1s3;
    std::array<double, 3> f2s3;
    std::array<double, 3> b3;
    std::array<double, 4> f1s4;
    std::array<double, 4> f2s4;
    std::array<double, 4> f3s4;
    std::array<double, 4> b4;
    for (size_t i = 0; i < MI_verts.size(); i++) {
        const auto& MI_vert = MI_verts[i];
        switch (MI_vert.simplex_size) {
        case 2: // on tet edge
        {
            auto vId1 = MI_vert.simplex_vert_indices[0];
            auto vId2 = MI_vert.simplex_vert_indices[1];
            auto fId1 = MI_vert.material_indices[0];
            auto fId2 = MI_vert.material_indices[1];
            auto f1 = funcVals(fId1, vId1) - funcVals(fId2, vId1);
            auto f2 = funcVals(fId1, vId2) - funcVals(fId2, vId2);
            //
            compute_barycentric_coords(f1, f2, b2);
            MI_pts[i][0] = b2[0] * pts[vId1][0] + b2[1] * pts[vId2][0];
            MI_pts[i][1] = b2[0] * pts[vId1][1] + b2[1] * pts[vId2][1];
            MI_pts[i][2] = b2[0] * pts[vId1][2] + b2[1] * pts[vId2][2];
            break;
        }
        case 3: // on tet face
        {
            auto vId1 = MI_vert.simplex_vert_indices[0];
            auto vId2 = MI_vert.simplex_vert_indices[1];
            auto vId3 = MI_vert.simplex_vert_indices[2];
            auto fId1 = MI_vert.material_indices[0];
            auto fId2 = MI_vert.material_indices[1];
            auto fId3 = MI_vert.material_indices[2];
            // f1 - f2
            f1s3[0] = funcVals(fId1, vId1) - funcVals(fId2, vId1);
            f1s3[1] = funcVals(fId1, vId2) - funcVals(fId2, vId2);
            f1s3[2] = funcVals(fId1, vId3) - funcVals(fId2, vId3);
            // f2 - f3
            f2s3[0] = funcVals(fId2, vId1) - funcVals(fId3, vId1);
            f2s3[1] = funcVals(fId2, vId2) - funcVals(fId3, vId2);
            f2s3[2] = funcVals(fId2, vId3) - funcVals(fId3, vId3);
            //
            compute_barycentric_coords(f1s3, f2s3, b3);
            MI_pts[i][0] = b3[0] * pts[vId1][0] + b3[1] * pts[vId2][0] + b3[2] * pts[vId3][0];
            MI_pts[i][1] = b3[0] * pts[vId1][1] + b3[1] * pts[vId2][1] + b3[2] * pts[vId3][1];
            MI_pts[i][2] = b3[0] * pts[vId1][2] + b3[1] * pts[vId2][2] + b3[2] * pts[vId3][2];
            break;
        }
        case 4: // in tet cell
        {
            auto vId1 = MI_vert.simplex_vert_indices[0];
            auto vId2 = MI_vert.simplex_vert_indices[1];
            auto vId3 = MI_vert.simplex_vert_indices[2];
            auto vId4 = MI_vert.simplex_vert_indices[3];
            auto fId1 = MI_vert.material_indices[0];
            auto fId2 = MI_vert.material_indices[1];
            auto fId3 = MI_vert.material_indices[2];
            auto fId4 = MI_vert.material_indices[3];
            // f1 - f2
            f1s4[0] = funcVals(fId1, vId1) - funcVals(fId2, vId1);
            f1s4[1] = funcVals(fId1, vId2) - funcVals(fId2, vId2);
            f1s4[2] = funcVals(fId1, vId3) - funcVals(fId2, vId3);
            f1s4[3] = funcVals(fId1, vId4) - funcVals(fId2, vId4);
            // f2 - f3
            f2s4[0] = funcVals(fId2, vId1) - funcVals(fId3, vId1);
            f2s4[1] = funcVals(fId2, vId2) - funcVals(fId3, vId2);
            f2s4[2] = funcVals(fId2, vId3) - funcVals(fId3, vId3);
            f2s4[3] = funcVals(fId2, vId4) - funcVals(fId3, vId4);
            // f3 - f4
            f3s4[0] = funcVals(fId3, vId1) - funcVals(fId4, vId1);
            f3s4[1] = funcVals(fId3, vId2) - funcVals(fId4, vId2);
            f3s4[2] = funcVals(fId3, vId3) - funcVals(fId4, vId3);
            f3s4[3] = funcVals(fId3, vId4) - funcVals(fId4, vId4);
            //
            compute_barycentric_coords(f1s4, f2s4, f3s4, b4);
            MI_pts[i][0] = b4[0] * pts[vId1][0] + b4[1] * pts[vId2][0] + b4[2] * pts[vId3][0] +
                            b4[3] * pts[vId4][0];
            MI_pts[i][1] = b4[0] * pts[vId1][1] + b4[1] * pts[vId2][1] + b4[2] * pts[vId3][1] +
                            b4[3] * pts[vId4][1];
            MI_pts[i][2] = b4[0] * pts[vId1][2] + b4[1] * pts[vId2][2] + b4[2] * pts[vId3][2] +
                            b4[3] * pts[vId4][2];
            break;
        }
        case 1: // on tet vertex
            MI_pts[i] = pts[MI_vert.simplex_vert_indices[0]];
            break;
        default: break;
        }
    }
}

void compute_face_order_in_one_tet_MI(const MaterialInterface<3>& tet_cut_result,
    const std::vector<PolygonFace>& MI_faces,
    const Edge& edge,
    std::vector<std::pair<size_t, int>>& ordered_faces)
{
    // in a single tetrahedron, num face pairs == num faces incident to the edge
    size_t num_incident_faces = edge.face_edge_indices.size();
    ordered_faces.resize(num_incident_faces);
    //
    std::vector<bool> is_incident_faces(tet_cut_result.faces.size(), false);
    // map: tet face id --> MI-face id
    std::vector<size_t> MI_face_Id_of_face(tet_cut_result.faces.size(), MaterialInterface<3>::None);
    for (const auto& fId_eId_pair : edge.face_edge_indices) {
        size_t MI_face_id = fId_eId_pair.first;
        size_t face_id = MI_faces[MI_face_id].tet_face_indices[0].second;
        is_incident_faces[face_id] = true;
        MI_face_Id_of_face[face_id] = MI_face_id;
    }

    // map: material label -> cell id
    size_t max_material_label = 0;
    for (const auto& cell : tet_cut_result.cells) {
        if (cell.material_label > max_material_label) {
            max_material_label = cell.material_label;
        }
    }
    std::vector<size_t> cell_of_material(max_material_label+1, MaterialInterface<3>::None);
    for (size_t i = 0; i < tet_cut_result.cells.size() ; ++i) {
        cell_of_material[tet_cut_result.cells[i].material_label] = i;
    }


    // travel around the edge
    size_t MI_face_id1 = edge.face_edge_indices[0].first;
    size_t face_id1 = MI_faces[MI_face_id1].tet_face_indices[0].second;
    int face_sign = 1;
    size_t cell_id = cell_of_material[tet_cut_result.faces[face_id1].positive_material_label];
    size_t face_id2 = MaterialInterface<3>::None;
    for (size_t i = 0; i < num_incident_faces; i++) {
        ordered_faces[i] = std::make_pair(MI_face_id1, face_sign);
        // find next face
        for (const auto& fId : tet_cut_result.cells[cell_id].faces) {
            if (is_incident_faces[fId] && fId != face_id1) {
                face_id2 = fId;
            }
        }
        // assert: face_id2 != None
        // get sign of face2 and find next cell
        //.....
        if (cell_of_material[tet_cut_result.faces[face_id2].positive_material_label] == cell_id) {
            cell_id = cell_of_material[tet_cut_result.faces[face_id2].negative_material_label];
            face_sign = -1;
        } else {
            cell_id = cell_of_material[tet_cut_result.faces[face_id2].positive_material_label];
            face_sign = 1;
        }
        // update face1 and clear face2
        face_id1 = face_id2;
        MI_face_id1 = MI_face_Id_of_face[face_id2];
        face_id2 = MaterialInterface<3>::None;
    }
}

void build_next_vert(const std::vector<std::array<double, 3>>& pts,
    const std::vector<std::array<size_t, 4>>& tets,
    std::vector<size_t>& next_vert)
{
    next_vert.resize(pts.size(), std::numeric_limits<size_t>::max());
    for (const auto& tet : tets) {
        // find the smallest vertex of tet
        size_t min_id = 0;
        for (size_t i = 1; i < 4; ++i) {
            if (point_xyz_less(pts[tet[i]], pts[tet[min_id]])) {
                min_id = i;
            }
        }
        size_t min_vId = tet[min_id];
        //
        for (size_t i = 0; i < 4; ++i) {
            if (i != min_id) {
                next_vert[tet[i]] = min_vId;
            }
        }
    }
}
void compute_edge_intersection_order_MI(const MaterialInterface<3>& tet_cut_result,
    size_t v,
    size_t u,
    std::vector<size_t>& vert_indices)
{
    const auto& vertices = tet_cut_result.vertices;
    const auto& faces = tet_cut_result.faces;
    //
    std::array<bool,4> edge_flag {true, true, true, true};
    edge_flag[v] = false;
    edge_flag[u] = false;

    // find vertices on edge v->u, and index of v and u in tet_cut_result.vertices
    size_t v_id, u_id;
    std::vector<bool> is_on_edge_vu(vertices.size(), false);
    size_t num_vert_in_edge_vu = 0;
    size_t flag_count;
    size_t other_plane;
    for (size_t i = 0; i < vertices.size(); ++i) {
        flag_count = 0;
        other_plane = MaterialInterface<3>::None;
        const auto& vert = vertices[i];
        // vert.size() == 4
        for (int j = 0; j < 4; ++j) {
            if (vert[j] < 4) {  // 0,1,2,3 are tet boundary
                if (edge_flag[vert[j]]) {
                    ++flag_count;
                } else {
                    other_plane = vert[j];
                }
            }
        }
        if (flag_count == 2) {
            is_on_edge_vu[i] = true;
            if (other_plane == u) {  // current vertex is v
                v_id = i;
            } else if (other_plane == v) {  // current vertex is u
                u_id = i;
            } else {  // current vertex is in interior of edge v->u
                ++num_vert_in_edge_vu;
            }
        }
    }
    if (num_vert_in_edge_vu == 0) { // no intersection points in edge v->u
        vert_indices.push_back(v_id);
        vert_indices.push_back(u_id);
        return;
    }

    // find all faces on triangle v->u->w, w is a third tet vertex
    size_t pId;  // plane index of v->u->w, pick pId != v and pId != u
    for (size_t i = 0; i < 4; ++i) {
        if (edge_flag[i]) {
            pId = i;
            break;
        }
    }
    std::vector<size_t> faces_on_tri_vuw;
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& face = faces[i];
        if (face.positive_material_label == pId) { // face is on tet boundary pId
            faces_on_tri_vuw.push_back(i);
        }
    }

    // build edge-face connectivity on triangle v->u->w
    // map: edge (v1, v2) --> incident faces (f1, f2), f2 is None if (v1,v2) is on triangle boundary
    absl::flat_hash_map<std::pair<size_t, size_t>, std::pair<size_t,size_t>> faces_of_edge;
    for (auto fId : faces_on_tri_vuw) {
        const auto& face = faces[fId];
        size_t num_vert = face.vertices.size();
        for (size_t i = 0; i < num_vert; ++i) {
            size_t i_next = (i + 1) % num_vert;
            size_t vi = face.vertices[i];
            size_t vi_next = face.vertices[i_next];
            // add fId to edge (vi, vi_next)
            auto iter_inserted = faces_of_edge.try_emplace(
                std::make_pair(vi, vi_next),
                std::make_pair(fId, MaterialInterface<3>::None));
            if (!iter_inserted.second) {  // inserted before
                iter_inserted.first->second.second = fId;
            }
            // add fId to edge (vi_next, vi)
            iter_inserted = faces_of_edge.try_emplace(
                std::make_pair(vi_next, vi),
                std::make_pair(fId, MaterialInterface<3>::None));
            if (!iter_inserted.second) {  // inserted before
                iter_inserted.first->second.second = fId;
            }
        }
    }

    // find the face on triangle v->u->w:
    // 1. has vertex v
    // 2. has an edge on v->u
    size_t f_start;
    for (auto fId : faces_on_tri_vuw) {
        const auto& face = faces[fId];
        bool find_v = false;
        size_t count = 0;
        for (auto vi : face.vertices) {
            if (vi == v_id) {
                find_v = true;
            }
            if (is_on_edge_vu[vi]) {
                ++count;
            }
        }
        if (find_v && count == 2) {
            f_start = fId;
            break;
        }
    }

    // trace edge v->u
    vert_indices.reserve(num_vert_in_edge_vu + 2);
    vert_indices.push_back(v_id);
    std::vector<bool> visited_face(faces.size(), false);
    size_t v_curr = v_id;
    size_t f_curr = f_start;
    std::pair<size_t, size_t> edge_prev;
    std::pair<size_t, size_t> edge_next;
    std::pair<size_t, size_t> edge_on_vu;
    while (v_curr != u_id) {
        // clear edge_on_vu
        edge_on_vu.first = MaterialInterface<3>::None;
        edge_on_vu.second = MaterialInterface<3>::None;
        //
        const auto& face = faces[f_curr];
        size_t num_vert = face.vertices.size();
        // visit all edges of face, find edge_prev, edge_next and edge_on_vu
        for (size_t i = 0; i < num_vert; ++i) {
            size_t i_next = (i+1) % num_vert;
            size_t vi = face.vertices[i];
            size_t vi_next = face.vertices[i_next];
            if (is_on_edge_vu[vi] && !is_on_edge_vu[vi_next]) {
                auto& two_faces = faces_of_edge[std::make_pair(vi, vi_next)];
                auto other_face = (two_faces.first == f_curr) ? two_faces.second : two_faces.first;
                if (vi == v_id || (other_face != MaterialInterface<3>::None && visited_face[other_face])) {
                    edge_prev.first = vi;
                    edge_prev.second = vi_next;
                } else {
                    edge_next.first = vi;
                    edge_next.second = vi_next;
                }
            } else if (is_on_edge_vu[vi_next] && !is_on_edge_vu[vi]) {
                auto& two_faces = faces_of_edge[std::make_pair(vi, vi_next)];
                auto other_face = (two_faces.first == f_curr) ? two_faces.second : two_faces.first;
                if (vi_next == v_id || (other_face != MaterialInterface<3>::None && visited_face[other_face])) {
                    edge_prev.first = vi;
                    edge_prev.second = vi_next;
                } else {
                    edge_next.first = vi;
                    edge_next.second = vi_next;
                }
            } else if (is_on_edge_vu[vi] && is_on_edge_vu[vi_next]) {
                edge_on_vu.first = vi;
                edge_on_vu.second = vi_next;
            }
        }
        //
        if (edge_on_vu.first == MaterialInterface<3>::None) {
            // no edge of the face is on v->u
            // keep v_curr, update f_curr
            visited_face[f_curr] = true;
            auto& two_faces = faces_of_edge[edge_next];
            f_curr = (two_faces.first == f_curr) ? two_faces.second : two_faces.first;
        } else {
            // there is an edge of the face on v->u
            // update v_curr to be the next vert on edge v->u
            v_curr = (edge_on_vu.first == v_curr) ? edge_on_vu.second : edge_on_vu.first;
            vert_indices.push_back(v_curr);
            // update f_curr
            visited_face[f_curr] = true;
            auto& two_faces = faces_of_edge[edge_next];
            f_curr = (two_faces.first == f_curr) ? two_faces.second : two_faces.first;
        }
    }
}

void compute_passing_face_pair_MI(const MaterialInterface<3>& tet_cut_result,
    size_t v1,
    size_t v2,
    std::pair<size_t, int>& face_orient1,
    std::pair<size_t, int>& face_orient2)
{
    // find a face incident to edge v1 -> v2
    const auto& faces = tet_cut_result.faces;
    size_t incident_face_id;
    bool found_incident_face = false;
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& face = faces[i];
        size_t num_vert = face.vertices.size();
        for (size_t j = 0; j < num_vert; ++j) {
            size_t vId1 = face.vertices[j];
            size_t vId2 = face.vertices[(j+1) % num_vert];
            if ((vId1==v1 && vId2==v2) || (vId1==v2 && vId2 == v1)) {
                incident_face_id = i;
                found_incident_face = true;
                break;
            }
        }
        if (found_incident_face) {
            break;
        }
    }
    // assert: found_incident_face == true

    // map: material label -> cell id
    size_t max_material_label = 0;
    for (const auto& cell : tet_cut_result.cells) {
        if (cell.material_label > max_material_label) {
            max_material_label = cell.material_label;
        }
    }
    std::vector<size_t> cell_of_material(max_material_label+1, MaterialInterface<3>::None);
    for (size_t i = 0; i < tet_cut_result.cells.size() ; ++i) {
        cell_of_material[tet_cut_result.cells[i].material_label] = i;
    }
    // faces[incident_face_id].positive_material_label will be a tet boundary
    size_t cell_id = cell_of_material[faces[incident_face_id].negative_material_label];
    const auto& cell = tet_cut_result.cells[cell_id];

    // find the two faces
    // 1. bounding the cell
    // 2. passing vertex v1 or v2
    for (auto fId : cell.faces) {
        const auto& face = faces[fId];
        bool found_v1 = false;
        bool found_v2 = false;
        for (auto vId : face.vertices) {
            if (vId == v1) {
                found_v1 = true;
            } else if (vId == v2) {
                found_v2 = true;
            }
        }
        if (found_v1 && !found_v2) {
            // the current face passes v1 but not v2
            face_orient1.first = fId;
            face_orient1.second = (cell_of_material[face.positive_material_label] == cell_id) ? 1 : -1;
        } else if (!found_v1 && found_v2) {
            // the current face passes v2 but not v1
            face_orient2.first = fId;
            face_orient2.second = (cell_of_material[face.positive_material_label] == cell_id) ? 1 : -1;
        }
    }
}

void compute_passing_face_MI(const MaterialInterface<3>& tet_cut_result,
    size_t v,
    size_t u,
    std::pair<size_t, int>& face_orient)
{
    // find a face incident to edge v -> u
    const auto& faces = tet_cut_result.faces;
    size_t incident_face_id;
    bool found_incident_face = false;
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& face = faces[i];
        size_t num_vert = face.vertices.size();
        for (size_t j = 0; j < num_vert; ++j) {
            size_t vId1 = face.vertices[j];
            size_t vId2 = face.vertices[(j+1) % num_vert];
            if ((vId1==v && vId2==u) || (vId1==u && vId2 == v)) {
                incident_face_id = i;
                found_incident_face = true;
                break;
            }
        }
        if (found_incident_face) {
            break;
        }
    }
    // assert: found_incident_face == true

    // map: material label -> cell id
    size_t max_material_label = 0;
    for (const auto& cell : tet_cut_result.cells) {
        if (cell.material_label > max_material_label) {
            max_material_label = cell.material_label;
        }
    }
    std::vector<size_t> cell_of_material(max_material_label+1, MaterialInterface<3>::None);
    for (size_t i = 0; i < tet_cut_result.cells.size() ; ++i) {
        cell_of_material[tet_cut_result.cells[i].material_label] = i;
    }
    // faces[incident_face_id].positive_material_label will be a tet boundary
    size_t cell_id = cell_of_material[faces[incident_face_id].negative_material_label];
    const auto& cell = tet_cut_result.cells[cell_id];

    // find the face
    // 1. bounding the cell
    // 2. passing vertex v
    for (auto fId : cell.faces) {
        const auto& face = faces[fId];
        bool found_v = false;
        bool found_u = false;
        for (auto vId : face.vertices) {
            if (vId == v) {
                found_v = true;
            } else if (vId == u) {
                found_u = true;
            }
        }
        if (found_v && !found_u) {
            // the current face passes v but not u
            face_orient.first = fId;
            face_orient.second = (cell_of_material[face.positive_material_label] == cell_id) ? 1 : -1;
            return;
        }
    }
}
