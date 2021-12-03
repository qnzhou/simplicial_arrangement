#pragma once

#include "common.h"

namespace simplicial_arrangement {

// Forward declarations.
template <typename Scalar, int DIM>
class MaterialInterfaceBuilder;

namespace internal {

void add_material(MaterialInterfaceBuilder<double, 2>& builder, size_t material_index);
void add_material(MaterialInterfaceBuilder<Int, 2>& builder, size_t material_index);

void add_material(MaterialInterfaceBuilder<double, 3>& builder, size_t material_index);
void add_material(MaterialInterfaceBuilder<Int, 3>& builder, size_t material_index);

} // namespace internal
} // namespace simplicitial_arrangement
