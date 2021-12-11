#ifndef LOOK_UP_TABLE_H
#define LOOK_UP_TABLE_H

#include <simplicial_arrangement/simplicial_arrangement.h>
#include <nlohmann/json.hpp>

namespace simplicial_arrangement {

// one function look up table
extern std::unique_ptr<std::vector<Arrangement<3>>> one_func_lookup_table;

// two function look up table
extern std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> to_check_edge_table;
extern std::unique_ptr<std::vector<std::vector<Arrangement<3>>>> two_func_lookup_table;

//
extern bool use_lookup_table;

bool load_lookup_table();

bool load_one_func_lookup_table(const std::string &filename);
bool load_to_check_edge_table(const std::string &filename);
bool load_two_func_lookup_table(const std::string &filename);

void load_arrangement(Arrangement<3> &arrangement, const nlohmann::json &data);

template <typename T>
void load_vector(std::vector<T> &vec, const nlohmann::json &data);


template <typename T>
void load_vector(std::vector<T> &vec, const nlohmann::json &data)
{
    vec.resize(data.size());
    for (size_t i = 0; i < vec.size(); i++) {
        vec[i] = data[i].get<T>();
    }
}

} // namespace simplicial_arrangement

#endif // !LOOK_UP_TABLE_H
