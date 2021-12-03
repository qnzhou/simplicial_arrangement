#pragma once

#include "BSPNode.h"
#include "DisjointSets.h"
#include "common.h"
#include "add_material.h"

#include <simplicial_arrangement/material_interface.h>

#include <absl/container/flat_hash_map.h>

#include <algorithm>
#include <array>
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
        const size_t num_materials = m_materials.size();
        if (num_materials < 1) {
            logger().error(
                "At least one material is needed to generate a valid material interface");
            throw std::runtime_error("Insufficient materials!");
        }
        m_unique_materials.init(num_materials);

        initialize_simplex_boundary();
        initialize_material_interface();

        for (size_t i = 1; i < num_materials; i++) {
            internal::add_material(*this, i + DIM + 1);
        }
    }

    const Material<Scalar, DIM>& get_material(size_t index) const
    {
        if (index < DIM + 1) {
            return m_simplex_boundary[index];
        } else {
            return m_materials[index - DIM - 1];
        }
    }

    const MaterialInterface<DIM>& get_material_interface() const {
        return m_material_interface;
    }

    MaterialInterface<DIM>& get_material_interface() {
        return m_material_interface;
    }

private:
    void initialize_simplex_boundary()
    {
        // We use psudo-material to represent simplex boundary.  They should not
        // be used in math computations directly.
        m_simplex_boundary.resize(DIM + 1);
        if constexpr (DIM == 2) {
            m_simplex_boundary[0] = {INFINITE, 0, 0};
            m_simplex_boundary[1] = {0, INFINITE, 0};
            m_simplex_boundary[2] = {0, 0, INFINITE};
        } else {
            m_simplex_boundary[0] = {INFINITE, 0, 0, 0};
            m_simplex_boundary[1] = {0, INFINITE, 0, 0};
            m_simplex_boundary[2] = {0, 0, INFINITE, 0};
            m_simplex_boundary[3] = {0, 0, 0, INFINITE};
        }
    }

    void initialize_material_interface()
    {
        // Initialize a valid material interface using the first material.
        m_material_interface.vertices.resize(DIM + 1);
        m_material_interface.faces.resize(DIM + 1);
        m_material_interface.cells.resize(1);
        m_material_interface.cells[0].material_label = DIM + 1;

        if constexpr (DIM == 2) {
            m_material_interface.vertices[0] = {1, 2, 3};
            m_material_interface.vertices[1] = {2, 0, 3};
            m_material_interface.vertices[2] = {0, 1, 3};

            m_material_interface.faces[0].vertices = {1, 2};
            m_material_interface.faces[1].vertices = {2, 0};
            m_material_interface.faces[2].vertices = {0, 1};

            m_material_interface.faces[0].positive_cell = 0;
            m_material_interface.faces[0].negative_cell = 3;
            m_material_interface.faces[1].positive_cell = 1;
            m_material_interface.faces[1].negative_cell = 3;
            m_material_interface.faces[2].positive_cell = 2;
            m_material_interface.faces[2].negative_cell = 3;
        } else {
            m_material_interface.vertices[0] = {1, 2, 3, 4};
            m_material_interface.vertices[1] = {0, 2, 3, 4};
            m_material_interface.vertices[2] = {0, 1, 3, 4};
            m_material_interface.vertices[2] = {0, 1, 2, 4};

            m_material_interface.faces[0].vertices = {1, 2, 3};
            m_material_interface.faces[1].vertices = {0, 3, 2};
            m_material_interface.faces[2].vertices = {0, 1, 3};
            m_material_interface.faces[3].vertices = {0, 2, 1};

            m_material_interface.faces[0].positive_cell = 0;
            m_material_interface.faces[0].negative_cell = 4;
            m_material_interface.faces[1].positive_cell = 1;
            m_material_interface.faces[1].negative_cell = 4;
            m_material_interface.faces[2].positive_cell = 2;
            m_material_interface.faces[2].negative_cell = 4;
            m_material_interface.faces[3].positive_cell = 3;
            m_material_interface.faces[3].negative_cell = 4;
        }
    }

    inline void sort(Joint<DIM>& p) const
    {
        if constexpr (DIM == 2) {
            if (p[0] > p[1]) std::swap(p[0], p[1]);
            if (p[0] > p[2]) std::swap(p[0], p[2]);
            if (p[1] > p[2]) std::swap(p[1], p[2]);
        } else if (DIM == 3) {
            if (p[0] > p[1]) std::swap(p[0], p[1]);
            if (p[0] > p[2]) std::swap(p[0], p[2]);
            if (p[0] > p[3]) std::swap(p[0], p[3]);
            if (p[1] > p[2]) std::swap(p[1], p[2]);
            if (p[1] > p[3]) std::swap(p[1], p[3]);
            if (p[2] > p[3]) std::swap(p[2], p[3]);
        }
    }

private:
    std::vector<Material<Scalar, DIM>> m_simplex_boundary;
    const std::vector<Material<Scalar, DIM>>& m_materials;
    utils::DisjointSets m_unique_materials;
    MaterialInterface<DIM> m_material_interface;
};

} // namespace simplicial_arrangement
