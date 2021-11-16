#include "SimplicialArrangementBuilder.h"
#include "MaterialInterfaceBuilder.h"

namespace simplicial_arrangement {

template <typename Scalar, int DIM>
Arrangement<DIM> compute_arrangement_impl(const std::vector<Plane<Scalar, DIM>>& planes)
{
    SimplicialArrangementBuilder<Scalar, DIM> builder(planes);
    return builder.extract_arrangement();
}

Arrangement<2> compute_arrangement(const std::vector<Plane<double, 2>>& planes)
{
    return compute_arrangement_impl<double, 2>(planes);
}

Arrangement<2> compute_arrangement(const std::vector<Plane<Int, 2>>& planes)
{
    return compute_arrangement_impl<Int, 2>(planes);
}

Arrangement<3> compute_arrangement(const std::vector<Plane<double, 3>>& planes)
{
    return compute_arrangement_impl<double, 3>(planes);
}

Arrangement<3> compute_arrangement(const std::vector<Plane<Int, 3>>& planes)
{
    return compute_arrangement_impl<Int, 3>(planes);
}

template <typename Scalar, int DIM>
Arrangement<DIM> compute_material_interface_impl(const std::vector<Plane<Scalar, DIM>>& planes)
{
    MaterialInterfaceBuilder<Scalar, DIM> builder(planes);
    return builder.extract_material_interface();
}

Arrangement<2> compute_material_interface(const std::vector<Plane<double, 2>>& planes)
{
    return compute_material_interface_impl<double, 2>(planes);
}

Arrangement<2> compute_material_interface(const std::vector<Plane<Int, 2>>& planes)
{
    return compute_material_interface_impl<Int, 2>(planes);
}

Arrangement<3> compute_material_interface(const std::vector<Plane<double, 3>>& planes)
{
    return compute_material_interface_impl<double, 3>(planes);
}

Arrangement<3> compute_material_interface(const std::vector<Plane<Int, 3>>& planes)
{
    return compute_material_interface_impl<Int, 3>(planes);
}

} // namespace simplicial_arrangement
