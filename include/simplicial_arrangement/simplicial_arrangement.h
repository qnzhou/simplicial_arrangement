#pragma once

#include <absl/numeric/int128.h>

#include <array>
#include <vector>

namespace simplicial_arrangement {

using Int = absl::int128;

/**
 * A plane is defined by the barycentric plane equation:
 *     f0 * b0 + f1 * b1 + f2 * b2 + f3 * b3 = 0 // For 3D
 *     f0 * b0 + f1 * b1 + f2 * b2 = 0           // For 2D
 * where the b's are the barycentric variables, and f's are the
 * plane equation coefficients.  We store the f's for each plane.
 */
template <typename Scalar, int DIM>
using Plane = std::array<Scalar, DIM + 1>;

/**
 * A point is represented as the intersection of planes.  We store the index
 * of the plane here.
 */
template <int DIM>
using Point = std::array<size_t, DIM>;

template <int DIM>
struct Arrangement;

Arrangement<2> compute_arrangement(const std::vector<Plane<double, 2>>& planes);
Arrangement<2> compute_arrangement(const std::vector<Plane<Int, 2>>& planes);
Arrangement<3> compute_arrangement(const std::vector<Plane<double, 3>>& planes);
Arrangement<3> compute_arrangement(const std::vector<Plane<Int, 3>>& planes);

/**
 * A self-contained data structure for 2D or 3D arrangement representation.
 */
template <int DIM>
struct Arrangement
{
    static_assert(DIM == 2 || DIM == 3, "Only 2D and 3D arrangements are supported.");
    static constexpr size_t None = std::numeric_limits<size_t>::max();

    /**
     * A face represents a (DIM-1)-dimensional polytope.  For 2D and 3D, its
     * orientation is uniquely defined by ordering of its boundary vertices.
     */
    struct Face
    {
        /**
         * An ordered list of boundary vertices.
         *
         * In 3D, the face is always oriented counterclockwise when viewed from
         * the positive side of the supporting plane.
         *
         * In 2D, the face (aka edge) is oriented such that the positive side of
         * the supporting plane (aka line) is on the right.
         */
        std::vector<size_t> vertices;

        /**
         * Plane index of the supporting plane.
         */
        size_t supporting_plane = None;

        /**
         * The cell index on the positive and negative side of this face.
         */
        size_t positive_cell = None;
        size_t negative_cell = None;
    };

    /**
     * A cell is a DIM-dimensional polytope.
     */
    struct Cell
    {
        /**
         * A set of boundary face indices in no particular order.
         */
        std::vector<size_t> faces;

        /**
         * The orientation of each boundary face with respect to this cell.
         * `face_orientations[i] == true` means this cell is on the positive
         * side of faces[i].
         */
        std::vector<bool> face_orientations;

        /**
         * The orientation with respect to all planes.
         * `plane_orientations[i] == true` means this cell is on the positive
         * side of plane i.
         *
         * TODO: This field contains an super set of the data stored in
         * `face_orientations`.  We probably should just keep one of them.
         */
        std::vector<bool> plane_orientations;
    };

    std::vector<Point<DIM>> vertices;
    std::vector<Face> faces;
    std::vector<Cell> cells;

    std::vector<size_t> unique_plane_indices;
    std::vector<std::vector<size_t>> unique_planes;
    std::vector<bool> unique_plane_orientations;
};

} // namespace simplicial_arrangement
