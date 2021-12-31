#pragma once

#include <simplicial_arrangement/material_interface.h>

#include <vector>

namespace simplicial_arrangement {

template <typename Scalar, int DIM>
class MaterialRepo
{
public:
    static_assert(DIM == 2 || DIM == 3, "Only 2D and 3D material interface are supported");
    static_assert(std::is_same<Scalar, double>::value || std::is_same<Scalar, Int>::value,
        "Only double and 128bit int are supported as Scalar.");
    static constexpr Scalar INFINITE = std::numeric_limits<Scalar>::max();

    MaterialRepo(const std::vector<Material<Scalar, DIM>>& materials)
        : m_materials(materials)
    {
        const size_t num_materials = m_materials.size();
        if (num_materials < 1) {
            logger().error(
                "At least one material is needed to generate a valid material interface");
            throw std::runtime_error("Insufficient materials!");
        }

        initialize_simplex_boundary();
    }

    const Material<Scalar, DIM>& get_material(size_t index) const
    {
        if (index < DIM + 1) {
            return m_simplex_boundary[index];
        } else {
            return m_materials[index - DIM - 1];
        }
    }

    size_t get_num_materials() const {
        return m_materials.size();
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

private:
    std::vector<Material<Scalar, DIM>> m_simplex_boundary;
    const std::vector<Material<Scalar, DIM>>& m_materials;
};

} // namespace simplicial_arrangement
