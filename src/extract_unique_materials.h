#pragma once

#include "MaterialRepo.h"

namespace simplicial_arrangement {

template <typename Scalar, int DIM>
auto extract_unique_materials(const MaterialRepo<Scalar, DIM>& repo)
    -> std::tuple<std::vector<size_t>, std::vector<std::vector<size_t>>>
{
    const size_t num_materials = repo.get_num_materials() + DIM + 1;

    // Special cases for 1 or 2 materials.
    if (num_materials == DIM + 2) {
        // One material.
        if constexpr (DIM == 2) {
            return {{0, 1, 2, 3}, {{0}, {1}, {2}, {3}}};
        } else {
            return {{0, 1, 2, 3, 4}, {{0}, {1}, {2}, {3}, {4}}};
        }
    } else if (num_materials == DIM + 3) {
        // Two materials.
        if constexpr (DIM == 2) {
            bool same_material = repo.get_material(3) == repo.get_material(4);
            if (same_material) {
                return {{0, 1, 2, 3, 3}, {{0}, {1}, {2}, {3, 4}}};
            } else {
                return {{0, 1, 2, 3, 4}, {{0}, {1}, {2}, {3}, {4}}};
            }
        } else {
            bool same_material = repo.get_material(4) == repo.get_material(5);
            if (same_material) {
                return {{0, 1, 2, 3, 4, 4}, {{0}, {1}, {2}, {3}, {4, 5}}};
            } else {
                return {{0, 1, 2, 3, 4, 5}, {{0}, {1}, {2}, {3}, {4}, {5}}};
            }
        }
    }

    std::vector<size_t> material_indices(num_materials);
    std::iota(material_indices.begin(), material_indices.end(), 0);

    auto material_less = [&](size_t i0, size_t i1) {
        const auto& m0 = repo.get_material(i0);
        const auto& m1 = repo.get_material(i1);
        for (size_t i = 0; i <= DIM; i++) {
            if (m0[i] < m1[i]) return true;
            if (m0[i] > m1[i]) return false;
        }
        return false;
    };

    std::sort(material_indices.begin(), material_indices.end(), material_less);

    std::vector<size_t> unique_material_indices;
    std::vector<std::vector<size_t>> unique_materials;
    unique_material_indices.resize(num_materials);
    size_t unique_id = 0;
    size_t id = 0;

    while (id < num_materials) {
        size_t i = material_indices[id];
        const auto& mi = repo.get_material(i);
        for (; id < num_materials; id++) {
            size_t j = material_indices[id];
            const auto& mj = repo.get_material(j);
            if (mi != mj) {
                break;
            }
            unique_material_indices[j] = unique_id;
        }
        unique_id++;
    }

    unique_materials.resize(unique_id);
    for (size_t i = 0; i < num_materials; i++) {
        size_t material_id = unique_material_indices[i];
        unique_materials[material_id].push_back(i);
    }

    logger().debug("Detected {} unique materials", unique_id);
    return {unique_material_indices, unique_materials};
}

} // namespace simplicial_arrangement
