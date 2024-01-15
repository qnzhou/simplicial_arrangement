#pragma once

#include <cassert>
#include <cstddef>
#include <numeric>
#include <vector>

namespace simplicial_arrangement::utils {

class DisjointSets
{
public:
    DisjointSets() = default;
    explicit DisjointSets(size_t n) { init(n); }

    void init(size_t n)
    {
        m_parent.resize(n);
        std::iota(m_parent.begin(), m_parent.end(), size_t(0));
    }

    size_t size() const { return m_parent.size(); }

    void clear() { m_parent.clear(); }

    size_t find(size_t i);

    size_t merge(size_t i, size_t j)
    {
        const size_t root_i = find(i);
        const size_t root_j = find(j);
        return m_parent[root_j] = root_i;
    }

    std::tuple<std::vector<std::vector<size_t>>, std::vector<size_t>> extract_disjoint_sets();
    size_t extract_disjoint_set_indices(std::vector<size_t>& index_map);

private:
    std::vector<size_t> m_parent;
};

} // namespace simplicial_arrangement::utils
