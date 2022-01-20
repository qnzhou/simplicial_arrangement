#include <simplicial_arrangement/material_interface.h>
#include "MaterialInterfaceBuilder.h"

namespace simplicial_arrangement {

template <typename Scalar, int DIM>
MaterialInterface<DIM> compute_material_interface_impl(const std::vector<Material<Scalar, DIM>>& materials)
{
    return MaterialInterfaceBuilder<Scalar, DIM>(materials).export_material_interface();
}

MaterialInterface<2> compute_material_interface(const std::vector<Material<double, 2>>& materials)
{
    return compute_material_interface_impl<double, 2>(materials);
}

#ifndef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
MaterialInterface<2> compute_material_interface(const std::vector<Material<Int, 2>>& materials)
{
    return compute_material_interface_impl<Int, 2>(materials);
}
#endif

MaterialInterface<3> compute_material_interface(const std::vector<Material<double, 3>>& materials)
{
    return compute_material_interface_impl<double, 3>(materials);
}

#ifndef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
MaterialInterface<3> compute_material_interface(const std::vector<Material<Int, 3>>& materials)
{
    return compute_material_interface_impl<Int, 3>(materials);
}
#endif

} // namespace simplicial_arrangement
