#pragma once

#include <vector>

namespace simplicial_arrangement {

/**
 * A self-contained data structure for 2D or 3D arrangement representation.
 */
template <int DIM>
struct Arrangement
{
    static_assert(DIM == 2 || DIM == 3, "Only 2D and 3D arrangements are supported.");
    static constexpr size_t None = std::numeric_limits<size_t>::max();

    /**
     * A vertex is represented by DIM intersecting hyperplanes.
     */
    using Vertex = std::array<size_t, DIM>;

    /**
     * A face represents a (DIM-1)-dimentional polytope.  For 2D and 3D, its
     * orientation is uniquely defined by ordering of its boundary vertices.
     */
    struct Face
    {
        /**
         * An ordered list of boundary vertices.  The ordering determines the
         * face orientation.
         */
        std::vector<size_t> vertices;

        /**
         * A set of supporing plane indices for this face.
         */
        std::vector<size_t> supporting_planes;

        /**
         * Whether this face is consistently orirented with each supporting
         * plane.
         */
        std::vector<bool> supporting_plane_orientations;

        /**
         * The cell index on the positive and negative side of this face.
         */
        size_t positive_cell = None;
        size_t negative_cell = None;
    };

    /**
     * A cell is a DIM-dimentional polytope.
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
         * side of face i.
         */
        std::vector<bool> face_orientations;
    };

    std::vector<Vertex> vertices;
    std::vector<Face> faces;
    std::vector<Cell> cells;
};

} // namespace simplicial_arrangement
