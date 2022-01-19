#pragma once

#include "common.h"
#include "robust_assert.h"

#include <implicit_predicates/implicit_predicates.h>

#include <vector>

namespace simplicial_arrangement::utils {

inline int8_t signof(implicit_predicates::Orientation o)
{
    ROBUST_ASSERT(o != implicit_predicates::INVALID);
    if (o > 0)
        return 1;
    else if (o < 0)
        return -1;
    else
        return 0;
}

template <typename T, typename KeepFunc>
std::vector<size_t> shrink(std::vector<T>& c, const KeepFunc& keep)
{
    const size_t s = c.size();
    std::vector<size_t> index_map(s, INVALID);
    size_t active_count = 0;
    for (size_t i = 0; i < s; i++) {
        if (!keep(i)) continue;

        if (i != active_count) {
            std::swap(c[active_count], c[i]);
        }
        index_map[i] = active_count;
        active_count++;
    }
    c.resize(active_count);

    return index_map;
}

template <typename EdgeType>
bool edges_are_ordered(const std::vector<EdgeType>& edges, const std::vector<size_t>& edge_ids)
{
    const size_t num_involved_edges = edge_ids.size();
    for (size_t i = 0; i < num_involved_edges; i++) {
        const auto& curr_e = edges[edge_ids[i]];
        const auto& next_e = edges[edge_ids[(i + 1) % num_involved_edges]];

        if (curr_e.vertices[0] == next_e.vertices[0] || curr_e.vertices[0] == next_e.vertices[1])
            continue;
        if (curr_e.vertices[1] == next_e.vertices[0] || curr_e.vertices[1] == next_e.vertices[1])
            continue;
        return false;
    }
    return true;
}

} // namespace simplicial_arrangement::utils
