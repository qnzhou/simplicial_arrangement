#pragma once

#include "common.h"

#include <vector>

namespace simplicial_arrangement {

template <typename T, typename KeepFunc>
std::vector<size_t> shrink(std::vector<T>& c, const KeepFunc& keep) {
    const size_t s = c.size();
    std::vector<size_t> index_map(s, INVALID);
    size_t active_count = 0;
    for (size_t i=0; i<s; i++) {
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

}
