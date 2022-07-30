#pragma once
#include <simplicial_arrangement/material_interface.h>
#include <simplicial_arrangement/simplicial_arrangement.h>

#include <memory>

namespace simplicial_arrangement {

// one function look up table
extern std::unique_ptr<std::vector<Arrangement<3>>> one_func_lookup_table;

// two function look up table
extern std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> to_check_edge_table;
extern std::unique_ptr<std::vector<std::vector<Arrangement<3>>>> two_func_lookup_table;

// Lookup table for simplicial arrangement
extern std::vector<Arrangement<3>> ar_data;
extern std::vector<size_t> ar_indices;

// Lookup table for material interface
extern std::vector<MaterialInterface<3>> mi_data;
extern std::vector<size_t> mi_indices;


extern bool use_lookup_table;


// For 1 plane arrangement lookup.
size_t ar_compute_outer_index(const Plane<Int, 3>& p0);
size_t ar_compute_outer_index(const Plane<double, 3>& p0);

// For 2 plane arrangement lookup.
size_t ar_compute_outer_index(const Plane<Int, 3>& p0, const Plane<Int, 3>& p1);
size_t ar_compute_outer_index(const Plane<double, 3>& p0, const Plane<double, 3>& p1);
size_t ar_compute_inner_index(
    size_t outer_index, const Plane<Int, 3>& p0, const Plane<Int, 3>& p1);
size_t ar_compute_inner_index(
    size_t outer_index, const Plane<double, 3>& p0, const Plane<double, 3>& p1);


// For 2 material lookup.
size_t mi_compute_outer_index(const Material<Int, 3>& m0, const Material<Int, 3>& m1);
size_t mi_compute_outer_index(const Material<double, 3>& m0, const Material<double, 3>& m1);

// For 3 material lookup.
size_t mi_compute_outer_index(
    const Material<Int, 3>& m0, const Material<Int, 3>& m1, const Material<Int, 3>& m2);
size_t mi_compute_outer_index(
    const Material<double, 3>& m0, const Material<double, 3>& m1, const Material<double, 3>& m2);
size_t mi_compute_inner_index(size_t outer_index,
    const Material<Int, 3>& m0,
    const Material<Int, 3>& m1,
    const Material<Int, 3>& m2);
size_t mi_compute_inner_index(size_t outer_index,
    const Material<double, 3>& m0,
    const Material<double, 3>& m1,
    const Material<double, 3>& m2);

} // namespace simplicial_arrangement
