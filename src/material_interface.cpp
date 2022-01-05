#include <simplicial_arrangement/material_interface.h>
#include "MaterialInterfaceBuilder.h"

namespace simplicial_arrangement {

template <typename Scalar, int DIM>
MaterialInterface<DIM> compute_material_interface_impl(const std::vector<Material<Scalar, DIM>>& materials)
{
    MaterialInterfaceBuilder<Scalar, DIM> builder(materials);
    return builder.get_material_interface();
}

MaterialInterface<2> compute_material_interface(const std::vector<Material<double, 2>>& materials)
{
    return compute_material_interface_impl<double, 2>(materials);
}

MaterialInterface<2> compute_material_interface(const std::vector<Material<Int, 2>>& materials)
{
    return compute_material_interface_impl<Int, 2>(materials);
}

MaterialInterface<3> compute_material_interface(const std::vector<Material<double, 3>>& materials)
{
    return compute_material_interface_impl<double, 3>(materials);
}

MaterialInterface<3> compute_material_interface(const std::vector<Material<Int, 3>>& materials)
{
    return compute_material_interface_impl<Int, 3>(materials);
}

} // namespace simplicial_arrangement