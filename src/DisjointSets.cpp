#include "DisjointSets.h"

namespace simplicial_arrangement::utils {

size_t DisjointSets::find(size_t i)
{
    while (m_parent[i] != i) {
        m_parent[i] = m_parent[m_parent[i]];
        i = m_parent[i];
    }
    return i;
}

std::tuple<std::vector<std::vector<size_t>>, std::vector<size_t>>
DisjointSets::extract_disjoint_sets()
{
    const size_t num_entries = size();
    std::vector<size_t> index_map;
    const size_t counter = extract_disjoint_set_indices(index_map);

    std::vector<std::vector<size_t>> disjoint_sets(counter);
    for (size_t i = 0; i < num_entries; i++) {
        disjoint_sets[index_map[i]].push_back(i);
    }
    return {disjoint_sets, index_map};
}

size_t DisjointSets::extract_disjoint_set_indices(std::vector<size_t>& index_map)
{
    constexpr size_t INVALID_VAL = std::numeric_limits<size_t>::max();
    const size_t num_entries = size();
    index_map.resize(num_entries, INVALID_VAL);
    size_t counter = 0;

    // Assign each roots a unique index.
    for (size_t i = 0; i < num_entries; i++) {
        const auto root = find(i);
        if (i == root) {
            index_map[i] = counter;
            counter++;
        }
    }

    // Assign all members the same index as their root.
    for (size_t i = 0; i < num_entries; i++) {
        const auto root = find(i);
        assert(index_map[root] != INVALID_VAL);
        index_map[i] = index_map[root];
    }

    return counter;
}

} // namespace simplicial_arrangement::utils
