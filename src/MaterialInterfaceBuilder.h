#pragma once

#include "BSPNode.h"
#include "DisjointSets.h"
#include "MIComplex.h"
#include "MaterialRepo.h"
#include "add_material.h"
#include "common.h"
#include "extract_material_interface.h"
#include "initialize_simplicial_mi_complex.h"

#include <simplicial_arrangement/material_interface.h>

#include <absl/container/flat_hash_map.h>

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
        const size_t num_materials = m_materials.get_num_materials();
        auto mi_complex = initialize_simplicial_mi_complex<DIM>();
        for (size_t i = 1; i < num_materials; i++) {
            internal::add_material(*this, mi_complex, i + DIM + 1);
        }
        m_material_interface = extract_material_interface(std::move(mi_complex));
    }

    const Material<Scalar, DIM>& get_material(size_t index) const
    {
        return m_materials.get_material(index);
    }

    const MaterialRepo<Scalar, DIM>& get_material_repo() const { return m_materials; }

    MaterialInterface<DIM>& get_material_interface() { return m_material_interface; }
    const MaterialInterface<DIM>& get_material_interface() const { return m_material_interface; }

private:
    MaterialRepo<Scalar, DIM> m_materials;
    utils::DisjointSets m_unique_materials;
    MaterialInterface<DIM> m_material_interface;
};

} // namespace simplicial_arrangement
