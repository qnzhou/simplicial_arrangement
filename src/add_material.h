#pragma once

#include "common.h"

namespace simplicial_arrangement {

// Forward declarations.
template <typename Scalar, int DIM>
class MaterialInterfaceBuilder;

template <int DIM>
struct MIComplex;

namespace internal {

void add_material(
    MaterialInterfaceBuilder<double, 2>& builder, MIComplex<2>& mi_complex, size_t material_index);
#ifndef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
void add_material(
    MaterialInterfaceBuilder<Int, 2>& builder, MIComplex<2>& mi_complex, size_t material_index);
#endif

void add_material(
    MaterialInterfaceBuilder<double, 3>& builder, MIComplex<3>& mi_complex, size_t material_index);
#ifndef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
void add_material(
    MaterialInterfaceBuilder<Int, 3>& builder, MIComplex<3>& mi_complex, size_t material_index);
#endif

} // namespace internal
} // namespace simplicial_arrangement
