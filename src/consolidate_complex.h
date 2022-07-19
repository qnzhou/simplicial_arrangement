#pragma once

#include <vector>
#include "utils.h"
#include "robust_assert.h"

namespace simplicial_arrangement {

/**
 * Remove unused verices and faces.
 */
template <template <int> typename Complex, int DIM>
void consolidate(Complex<DIM>& data)
{
    // Shrink faces.
    if constexpr (DIM == 3) {
        std::vector<bool> active_faces(data.faces.size(), false);
        for (auto& c : data.cells) {
            for (auto fid : c.faces) {
                active_faces[fid] = true;
            }
        }

        auto index_map = utils::shrink(data.faces, [&](size_t fid) { return active_faces[fid]; });

        for (auto& c : data.cells) {
            std::transform(c.faces.begin(), c.faces.end(), c.faces.begin(), [&](size_t i) {
                ROBUST_ASSERT(index_map[i] != INVALID);
                return index_map[i];
            });
        }
    }

    // Shrink edges.
    {
        std::vector<bool> active_edges(data.edges.size(), false);
        for (auto& f : data.faces) {
            for (auto eid : f.edges) {
                active_edges[eid] = true;
            }
        }

        auto index_map = utils::shrink(data.edges, [&](size_t eid) { return active_edges[eid]; });

        for (auto& f : data.faces) {
            std::transform(f.edges.begin(), f.edges.end(), f.edges.begin(), [&](size_t i) {
                ROBUST_ASSERT(index_map[i] != INVALID);
                return index_map[i];
            });
        }
    }

    // Shrink vertices.
    {
        std::vector<bool> active_vertices(data.vertices.size(), false);
        for (auto& e : data.edges) {
            for (auto vid : e.vertices) {
                ROBUST_ASSERT(vid != INVALID);
                active_vertices[vid] = true;
            }
        }

        auto index_map =
            utils::shrink(data.vertices, [&](size_t vid) { return active_vertices[vid]; });

        for (auto& e : data.edges) {
            e.vertices[0] = index_map[e.vertices[0]];
            e.vertices[1] = index_map[e.vertices[1]];
        }
    }
}

}
