#pragma once
#include <simplicial_arrangement/material_interface.h>
#include <simplicial_arrangement/simplicial_arrangement.h>

namespace simplicial_arrangement {

// one function look up table
extern std::unique_ptr<std::vector<Arrangement<3>>> one_func_lookup_table;

// two function look up table
extern std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> to_check_edge_table;
extern std::unique_ptr<std::vector<std::vector<Arrangement<3>>>> two_func_lookup_table;

// Lookup table for material interface
extern std::vector<MaterialInterface<3>> mi_data;
extern std::vector<size_t> mi_indices;

extern bool use_lookup_table;

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
