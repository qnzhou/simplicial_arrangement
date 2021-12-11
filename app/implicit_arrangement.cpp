#include <simplicial_arrangement/simplicial_arrangement.h>
#include <simplicial_arrangement/lookup_table.h>

#include <stdlib.h>
#include <iostream>

using namespace simplicial_arrangement;

bool load_tet_mesh(const std::string& filename,
    std::vector<std::vector<double>>& pts,
    std::vector<std::vector<size_t>>& tets)
{
    //todo
    return true;
}

double sphere_function(const std::vector<double> center, double r, const std::vector<double> p) {
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

    // load tet mesh
    std::string tet_mesh_file = "";
    std::vector<std::vector<double>> pts;
    std::vector<std::vector<size_t>> tets;
    load_tet_mesh(tet_mesh_file, pts, tets);

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
    std::vector<Arrangement<3>> cut_results(tets.size());
    for (size_t i = 0; i < tets.size(); i++) {
        const auto& func_ids = func_in_tet[i];
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

    // extract arrangement mesh
    // first try: without removing duplicate vertices and faces



    return 0;
}
