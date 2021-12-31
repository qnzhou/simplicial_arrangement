#pragma once

#include "BSPNode.h"
#include "DisjointSets.h"
#include "MIComplex.h"
#include "MaterialRepo.h"
#include "add_material.h"
#include "common.h"
#include "extract_material_interface.h"

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
        const size_t num_materials = m_materials.get_num_materials();
        auto mi_complex = initialize_complex();
        for (size_t i = 1; i < num_materials; i++) {
            internal::add_material(*this, mi_complex, i + DIM + 1);
        }
        m_material_interface = extract_material_interface(std::move(mi_complex));
    }

    const Material<Scalar, DIM>& get_material(size_t index) const
    {
        return m_materials.get_material(index);
    }

    const MaterialRepo<Scalar, DIM>& get_material_repo() const
    {
        return m_materials;
    }

    MaterialInterface<DIM>& get_material_interface() { return m_material_interface; }
    const MaterialInterface<DIM>& get_material_interface() const { return m_material_interface; }

private:
    MIComplex<DIM> initialize_complex()
    {
        MIComplex<DIM> mi_complex;

        // Initialize a valid material interface using the first material.
        mi_complex.vertices.resize(DIM + 1);

        if constexpr (DIM == 2) {
            mi_complex.vertices[0] = {1, 2, 3};
            mi_complex.vertices[1] = {2, 0, 3};
            mi_complex.vertices[2] = {0, 1, 3};

            mi_complex.edges.resize(3);
            mi_complex.edges[0].vertices = {1, 2};
            mi_complex.edges[1].vertices = {2, 0};
            mi_complex.edges[2].vertices = {0, 1};

            mi_complex.edges[0].positive_material_label = 0;
            mi_complex.edges[0].negative_material_label = 3;
            mi_complex.edges[1].positive_material_label = 1;
            mi_complex.edges[1].negative_material_label = 3;
            mi_complex.edges[2].positive_material_label = 2;
            mi_complex.edges[2].negative_material_label = 3;

            mi_complex.faces.resize(1);
            mi_complex.faces[0].material_label = 3;
            mi_complex.faces[0].edges = {0, 1, 2};
        } else {
            mi_complex.vertices[0] = {1, 2, 3, 4};
            mi_complex.vertices[1] = {0, 2, 3, 4};
            mi_complex.vertices[2] = {0, 1, 3, 4};
            mi_complex.vertices[3] = {0, 1, 2, 4};

            mi_complex.edges.resize(6);
            mi_complex.edges[0].vertices = {0, 1};
            mi_complex.edges[1].vertices = {0, 2};
            mi_complex.edges[2].vertices = {0, 3};
            mi_complex.edges[3].vertices = {1, 2};
            mi_complex.edges[4].vertices = {1, 3};
            mi_complex.edges[5].vertices = {2, 3};

            mi_complex.edges[0].supporting_materials = {2, 3, 4};
            mi_complex.edges[1].supporting_materials = {1, 3, 4};
            mi_complex.edges[2].supporting_materials = {1, 2, 4};
            mi_complex.edges[3].supporting_materials = {0, 3, 4};
            mi_complex.edges[4].supporting_materials = {0, 2, 4};
            mi_complex.edges[5].supporting_materials = {0, 1, 4};

            mi_complex.faces.resize(4);
            mi_complex.faces[0].edges = {3, 5, 4};
            mi_complex.faces[1].edges = {1, 2, 5};
            mi_complex.faces[2].edges = {0, 4, 2};
            mi_complex.faces[3].edges = {0, 1, 3};

            mi_complex.faces[0].positive_material_label = 0;
            mi_complex.faces[0].negative_material_label = 4;
            mi_complex.faces[1].positive_material_label = 1;
            mi_complex.faces[1].negative_material_label = 4;
            mi_complex.faces[2].positive_material_label = 2;
            mi_complex.faces[2].negative_material_label = 4;
            mi_complex.faces[3].positive_material_label = 3;
            mi_complex.faces[3].negative_material_label = 4;

            mi_complex.cells.resize(1);
            mi_complex.cells[0].material_label = 4;
            mi_complex.cells[0].faces = {0, 1, 2, 3};
        }
        return mi_complex;
    }

private:
    MaterialRepo<Scalar, DIM> m_materials;
    utils::DisjointSets m_unique_materials;
    MaterialInterface<DIM> m_material_interface;
};

} // namespace simplicial_arrangement
