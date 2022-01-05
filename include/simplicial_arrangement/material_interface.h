#pragma once
#include <absl/numeric/int128.h>

#include <array>
#include <vector>

namespace simplicial_arrangement {

using Int = absl::int128;

/**
 * A material is defined by the barycentric function:
 *     f(b0,b1,b2,b3) = m0 * b0 + m1 * b1 + m2 * b2 + m3 * b3   // For 3D
 *        f(b0,b1,b2) = m0 * b0 + m1 * b1 + m2 * b2             // For 2D
 * where the b's are the barycentric variables, and m's are the material value
 * at the vertices of a simplex.  We store the m's for each material.
 */
template <typename Scalar, int DIM>
using Material = std::array<Scalar, DIM + 1>;

/**
 * A material interface joint is where DIM+1 materials share the same value:
 *   f_0(p) = f_1(p) = ... = f_DIM(p)
 * We store the index of involved material labels here.
 */
template <int DIM>
using Joint = std::array<size_t, DIM + 1>;

template <int DIM>
struct MaterialInterface;

MaterialInterface<2> compute_material_interface(const std::vector<Material<double, 2>>& materials);
MaterialInterface<2> compute_material_interface(const std::vector<Material<Int, 2>>& materials);
MaterialInterface<3> compute_material_interface(const std::vector<Material<double, 3>>& materials);
MaterialInterface<3> compute_material_interface(const std::vector<Material<Int, 3>>& materials);

/**
 * A self-contained data structure for 2D and 3D material interface
 * representation.
 */
template <int DIM>
struct MaterialInterface
{
    static_assert(DIM == 2 || DIM == 3, "Only 2D and 3D material interface are supported.");
    static constexpr size_t None = std::numeric_limits<size_t>::max();

    /**
     * A face represents a (DIM-1)-dimensional polytope.  For 2D and 3D, its
     * orientation is uniquely defined by the ordering of its boundary vertices.
     */
    struct Face
    {
        std::vector<size_t> vertices; // ordered.
        size_t positive_material_label = None;
        size_t negative_material_label = None;
    };

    /**
     * A cell is a DIM-dimensional polytope.
     */
    struct Cell
    {
        std::vector<size_t> faces;
        size_t material_label = None;
    };

    std::vector<Joint<DIM>> vertices;
    std::vector<Face> faces;
    std::vector<Cell> cells;

    std::vector<size_t> unique_material_indices;
    std::vector<std::vector<size_t>> unique_materials;
};

} // namespace simplicial_arrangement
