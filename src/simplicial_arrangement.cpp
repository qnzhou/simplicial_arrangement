#include <implicit_predicates/implicit_predicates.h>
#include <simplicial_arrangement/lookup_table.h>

#include "ArrangementBuilder.h"
#include "SimplicialArrangementBuilder.h"
#include "MaterialInterfaceBuilder.h"
#include "lookup_table_forward_declarations.h"

namespace simplicial_arrangement {

template <typename Scalar, int DIM>
Arrangement<DIM> compute_arrangement_impl(const std::vector<Plane<Scalar, DIM>>& planes)
{
    //SimplicialArrangementBuilder<Scalar, DIM> builder(planes);
    //return builder.extract_arrangement();

    return ArrangementBuilder<Scalar, DIM>(planes).export_arrangement();
}

template <typename Scalar>
Arrangement<3> compute_arrangement_1_plane(const std::vector<Plane<Scalar, 3>>& planes)
{
    const Plane<Scalar, 3>& plane = planes[0];
    bool has_zero = false;
    size_t index = 0;
    for (size_t i = 0; i < 4; i++) {
        index <<= 1;
        if (plane[i] > 0) {
            index += 1;
        } else if (plane[i] == 0) {
            has_zero = true;
            break;
        }
    }

    if (has_zero) {
        return compute_arrangement_impl<Scalar, 3>(planes);
    } else { // return index-th entry of the look up table
        return (*one_func_lookup_table)[index];
    }
}

template <typename Scalar>
Arrangement<3> compute_arrangement_2_planes(const std::vector<Plane<Scalar, 3>>& planes)
{
    const Plane<Scalar, 3>& plane1 = planes[0];
    const Plane<Scalar, 3>& plane2 = planes[1];

    bool has_zero = false;
    size_t index1 = 0;
    for (size_t i = 0; i < 4; i++) {
        index1 <<= 1;
        if (plane1[i] > 0) {
            index1 += 1;
        } else if (plane1[i] == 0) {
            has_zero = true;
            break;
        }
    }
    if (has_zero) {
        return compute_arrangement_impl<Scalar, 3>(planes);
    }

    size_t index2 = 0;
    for (size_t i = 0; i < 4; i++) {
        index2 <<= 1;
        if (plane2[i] > 0) {
            index2 += 1;
        } else if (plane2[i] == 0) {
            has_zero = true;
            break;
        }
    }
    if (has_zero) {
        return compute_arrangement_impl<Scalar, 3>(planes);
    }

    //
    size_t index = (index1 << 4) + index2;
    // get orientation on edges
    const std::vector<std::pair<int, int>>& test_edges = (*to_check_edge_table)[index];
    size_t e_index = 0;
    for (const auto& edge : test_edges) {
        e_index <<= 1;
        Scalar vals1[2] = {plane1[edge.first], plane1[edge.second]};
        Scalar vals2[2] = {plane2[edge.first], plane2[edge.second]};
        auto orient = implicit_predicates::orient1d(vals1, vals2);
        if (orient == implicit_predicates::POSITIVE) {
            e_index += 1;
        } else if (orient == implicit_predicates::ZERO) {
            has_zero = true;
            break;
        }
    }
    if (has_zero) {
        return compute_arrangement_impl<Scalar, 3>(planes);
    }

    // take result from look-up table   [index][e_index]
    return (*two_func_lookup_table)[index][e_index];
}


template <typename Scalar>
Arrangement<3> compute_arrangement_lookup(const std::vector<Plane<Scalar, 3>>& planes)
{
    switch (planes.size()) {
    case 1: return compute_arrangement_1_plane<Scalar>(planes); break;
    case 2: return compute_arrangement_2_planes<Scalar>(planes); break;
    default: return compute_arrangement_impl<Scalar, 3>(planes);
    }
}

Arrangement<2> compute_arrangement(const std::vector<Plane<double, 2>>& planes)
{
    return compute_arrangement_impl<double, 2>(planes);
}

#ifndef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
Arrangement<2> compute_arrangement(const std::vector<Plane<Int, 2>>& planes)
{
    return compute_arrangement_impl<Int, 2>(planes);
}
#endif

Arrangement<3> compute_arrangement(const std::vector<Plane<double, 3>>& planes)
{
    if (use_lookup_table && two_func_lookup_table) {
        return compute_arrangement_lookup<double>(planes);
    }
    return compute_arrangement_impl<double, 3>(planes);
}

#ifndef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
Arrangement<3> compute_arrangement(const std::vector<Plane<Int, 3>>& planes)
{
    if (use_lookup_table && two_func_lookup_table) {
        return compute_arrangement_lookup<Int>(planes);
    }
    return compute_arrangement_impl<Int, 3>(planes);
}
#endif


} // namespace simplicial_arrangement
