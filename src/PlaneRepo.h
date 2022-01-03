#pragma once

#include "common.h"

#include <simplicial_arrangement/simplicial_arrangement.h>

#include <vector>

namespace simplicial_arrangement {

template <typename Scalar, int DIM>
class PlaneRepo
{
public:
    static_assert(DIM == 2 || DIM == 3, "Only 2D and 3D arrangement are supported");
    static_assert(std::is_same<Scalar, double>::value || std::is_same<Scalar, Int>::value,
        "Only double and 128bit int are supported as Scalar.");

    PlaneRepo(const std::vector<Plane<Scalar, DIM>>& planes)
        : m_planes(planes)
    {
        initialize_simplicial_planes();
    }

    const Plane<Scalar, DIM>& get_plane(size_t i) const
    {
        if (i <= DIM) {
            return m_simplex_planes[i];
        } else {
            return m_planes[i - DIM - 1];
        }
    }

    size_t get_num_planes() const { return m_planes.size(); }

private:
    void initialize_simplicial_planes()
    {
        if constexpr (DIM == 2) {
            m_simplex_planes[0] = {1, 0, 0};
            m_simplex_planes[1] = {0, 1, 0};
            m_simplex_planes[2] = {0, 0, 1};
        } else {
            m_simplex_planes[0] = {1, 0, 0, 0};
            m_simplex_planes[1] = {0, 1, 0, 0};
            m_simplex_planes[2] = {0, 0, 1, 0};
            m_simplex_planes[3] = {0, 0, 0, 1};
        }
    }

private:
    std::array<Plane<Scalar, DIM>, DIM + 1> m_simplex_planes;
    const std::vector<Plane<Scalar, DIM>>& m_planes;
};

} // namespace simplicial_arrangement
