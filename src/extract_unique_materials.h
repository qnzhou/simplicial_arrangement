#pragma once

#include "MaterialRepo.h"

namespace simplicial_arrangement {

template <typename Scalar, int DIM>
auto extract_unique_materials(const MaterialRepo<Scalar, DIM>& repo)
    -> std::tuple<std::vector<size_t>, std::vector<std::vector<size_t>>>
{
    const size_t num_materials = repo.get_num_materials() + DIM + 1;
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
