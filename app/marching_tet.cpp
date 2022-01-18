#include <simplicial_arrangement/simplicial_arrangement.h>
#include <simplicial_arrangement/lookup_table.h>

#include "implicit_arrangement_util.h"

#include <iostream>
#include "ScopedTimer.h"

using namespace simplicial_arrangement;

int main(int argc, const char* argv[]) {
    std::cout << "load table ..." << std::endl;
    bool loaded = load_lookup_table();
    if (loaded) {
        std::cout << "loading complete." << std::endl;
    }

    // load tet mesh
//    std::string dataDir = "D:/research/simplicial_arrangement/data/";
    std::string dataDir = "/Users/charlesdu/Downloads/research/implicit_modeling/code/simplicial_arrangement/data/";
    std::string resolution = "80k";
    std::string tet_mesh_file = dataDir + "tet_mesh_" + resolution + ".json";
    std::vector<std::array<double, 3>> pts;
    std::vector<std::array<size_t, 4>> tets;
    load_tet_mesh(tet_mesh_file, pts, tets);
    std::cout << "tet mesh: " << pts.size() << " verts, " << tets.size() << " tets." << std::endl;

    // load implicit function values, or evaluate
    size_t n_func = 4;
    std::vector<std::array<double, 3>> centers(n_func);
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
    std::vector<std::vector<bool>> has_isosurface(n_func);
    {
        ScopedTimer<> timer("filter out intersecting implicits in each tet");
        for (size_t i = 0; i < n_func; i++) {
            has_isosurface[i].resize(tets.size());
            for (size_t j = 0; j < tets.size(); j++) {
                int pos_count = 0;
                int neg_count = 0;
                for (const auto& vId : tets[j]) {
                    if (funcSigns[i][vId] == 1) {
                        pos_count += 1;
                    } else if (funcSigns[i][vId] == -1) {
                        neg_count += 1;
                    }
                }
                // tets[i].size() == 4
                if (pos_count < 4 && neg_count < 4) {
                    has_isosurface[i][j] = true; // function i is in tetrahedron j
                }                
            }
        }        
    }

    // marching tet on each implicit function
    std::vector<std::vector<Arrangement<3>>> cut_results(n_func);
    {
        ScopedTimer<> timer("marching tet");
        for (size_t i = 0; i < n_func; i++) {
            auto& marching_tet_results = cut_results[i];
            marching_tet_results.resize(tets.size());
            const auto& cur_has_isosurface = has_isosurface[i];
            const auto& cur_funcVals = funcVals[i];
            for (size_t j = 0; j < tets.size(); j++) {                
                if (cur_has_isosurface[j]) {
                    const auto& tet = tets[j];
                    std::vector<Plane<double, 3>> planes(1);
                    planes[0] = {cur_funcVals[tet[0]],
                        cur_funcVals[tet[1]],
                        cur_funcVals[tet[2]],
                        cur_funcVals[tet[3]]};
                    marching_tet_results[j] = compute_arrangement(planes);
                }
            }
        }
    }

    // extract arrangement mesh
    std::vector<std::vector<IsoVert>> iso_verts_list(n_func);
    std::vector<std::vector<PolygonFace>> iso_faces_list(n_func);
    {
        ScopedTimer<> timer("extract iso mesh (topology only)");
        for (size_t i = 0; i < n_func; i++) {
            extract_iso_mesh_marching_tet(
                has_isosurface[i], cut_results[i], tets, iso_verts_list[i], iso_faces_list[i]);
        }
    }

    // compute xyz of iso-vertices
    std::vector<std::vector<std::array<double, 3>>> iso_pts_list(n_func);
    {
        ScopedTimer<> timer("compute xyz of iso-vertices");
        for (size_t i = 0; i < n_func; i++) {
            compute_iso_vert_xyz_marching_tet(iso_verts_list[i], funcVals[i], pts, iso_pts_list[i]);
        }        
    }

    // test: export iso-mesh
    save_iso_mesh_list(dataDir + "marching_tet_iso_mesh_" + resolution + ".json",
        iso_pts_list,
        iso_faces_list);

    return 0;
}