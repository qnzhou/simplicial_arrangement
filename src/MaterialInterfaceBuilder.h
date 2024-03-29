#pragma once

#include "BSPNode.h"
#include "DisjointSets.h"
#include "MIComplex.h"
#include "MaterialRepo.h"
#include "add_material.h"
#include "common.h"
#include "extract_material_interface.h"
#include "extract_unique_materials.h"
#include "initialize_simplicial_mi_complex.h"
#include "lookup_table_forward_declarations.h"

#include <simplicial_arrangement/material_interface.h>

#include <absl/container/flat_hash_map.h>

#include <algorithm>
#include <type_traits>

namespace simplicial_arrangement {

template <typename Scalar, int DIM>
class MaterialInterfaceBuilder
{
public:
    static_assert(DIM == 2 || DIM == 3, "Only 2D and 3D material interface are supported");
    static_assert(std::is_same<Scalar, double>::value || std::is_same<Scalar, Int>::value,
        "Only double and 128bit int are supported as Scalar.");
    static constexpr Scalar INFINITE = std::numeric_limits<Scalar>::max();

public:
    MaterialInterfaceBuilder(const std::vector<Material<Scalar, DIM>>& materials)
        : m_materials(materials)
    {
        const size_t num_materials = materials.size();
        if (num_materials == 1) {
            auto mi_complex = initialize_simplicial_mi_complex<DIM>();
            m_material_interface = extract_material_interface(std::move(mi_complex));
            return;
        }

        bool need_full_compute = true;
        if (use_lookup_table) {
            const auto* mi = lookup(materials);
            if (mi != nullptr) {
                m_material_interface = *mi;
                need_full_compute = false;
            }
        }

        if (need_full_compute) {
            auto [unique_material_indices, unique_materials] =
                extract_unique_materials(m_materials);

            auto mi_complex = initialize_simplicial_mi_complex<DIM>();
            std::vector<bool> processed(unique_materials.size(), false);
            // material DIM+1 is already added by default.
            processed[unique_material_indices[DIM + 1]] = true;

            for (size_t i = 1; i < num_materials; i++) {
                const size_t material_index = i + DIM + 1;

                if (processed[unique_material_indices[material_index]]) continue;
                processed[unique_material_indices[material_index]] = true;

                internal::add_material(*this, mi_complex, material_index);
            }

            m_material_interface = extract_material_interface(std::move(mi_complex));
            if (unique_materials.size() != DIM + 1 + num_materials) {
                // Input materials contain duplicates.
                m_material_interface.unique_material_indices = std::move(unique_material_indices);
                m_material_interface.unique_materials = std::move(unique_materials);
            }
        }
    }

    const Material<Scalar, DIM>& get_material(size_t index) const
    {
        return m_materials.get_material(index);
    }

    const MaterialRepo<Scalar, DIM>& get_material_repo() const { return m_materials; }

    MaterialInterface<DIM>& get_material_interface() { return m_material_interface; }
    const MaterialInterface<DIM>& get_material_interface() const { return m_material_interface; }
    MaterialInterface<DIM>&& export_material_interface() && { return std::move(m_material_interface); }

private:
    const MaterialInterface<DIM>* lookup(const std::vector<Material<Scalar, DIM>>& materials) const
    {
        extern std::vector<MaterialInterface<3>> mi_data;
        extern std::vector<size_t> mi_indices;
        const size_t num_materials = materials.size();

        if constexpr (DIM == 3) {
            if (num_materials == 2) {
                const size_t outer_index = mi_compute_outer_index(materials[0], materials[1]);
                if (outer_index == INVALID) return nullptr;

                logger().debug("MI lookup outer index: {}", outer_index);
                const size_t start_idx = mi_indices[outer_index];
                assert(mi_indices[outer_index + 1] == start_idx + 1);

                logger().debug("MI lookup data index: {}", start_idx);
                return &mi_data[start_idx];
            } else if (num_materials == 3) {
                const size_t outer_index =
                    mi_compute_outer_index(materials[0], materials[1], materials[2]);
                if (outer_index == INVALID) return nullptr;

                logger().debug("MI lookup outer index: {}", outer_index);
                const size_t start_idx = mi_indices[outer_index];
                const size_t end_idx = mi_indices[outer_index + 1];
                if (end_idx == start_idx + 1) {
                    logger().debug("MI lookup data index: {}", start_idx);
                    return &mi_data[start_idx];
                } else if (end_idx > start_idx) {
                    const size_t inner_index = mi_compute_inner_index(
                        outer_index, materials[0], materials[1], materials[2]);
                    if (inner_index == INVALID) return nullptr;
                    assert(inner_index < end_idx - start_idx);
                    logger().debug("MI lookup inner index: {}", inner_index);
                    logger().debug("MI lookup data index: {}", start_idx + inner_index);
                    return &mi_data[start_idx + inner_index];
                }
            }
        }
        return nullptr;
    }

private:
    MaterialRepo<Scalar, DIM> m_materials;
    MaterialInterface<DIM> m_material_interface;
};

} // namespace simplicial_arrangement
