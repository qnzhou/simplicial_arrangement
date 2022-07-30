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
    return compute_arrangement_impl<double, 3>(planes);
}

#ifndef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
Arrangement<3> compute_arrangement(const std::vector<Plane<Int, 3>>& planes)
{
    return compute_arrangement_impl<Int, 3>(planes);
}
#endif


} // namespace simplicial_arrangement
