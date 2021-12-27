#pragma once

namespace simplicial_arrangement {

template <int DIM>
struct MaterialInterface;

template <int DIM>
struct MIComplex;

MaterialInterface<2> extract_material_interface(MIComplex<2>&& mi_complex);
MaterialInterface<3> extract_material_interface(MIComplex<3>&& mi_complex);

}
