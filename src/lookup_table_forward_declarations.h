#pragma once
#include <simplicial_arrangement/simplicial_arrangement.h>

namespace simplicial_arrangement {

// one function look up table
extern std::unique_ptr<std::vector<Arrangement<3>>> one_func_lookup_table;

// two function look up table
extern std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> to_check_edge_table;
extern std::unique_ptr<std::vector<std::vector<Arrangement<3>>>> two_func_lookup_table;

extern bool use_lookup_table;

} // namespace simplicial_arrangement
