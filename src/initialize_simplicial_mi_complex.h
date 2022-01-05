#pragma once

#include "MIComplex.h"

namespace simplicial_arrangement {

template <int DIM>
MIComplex<DIM> initialize_simplicial_mi_complex()
{
    MIComplex<DIM> mi_complex;

    // Initialize a valid material interface using the first material.
    mi_complex.vertices.resize(DIM + 1);

    if constexpr (DIM == 2) {
        mi_complex.vertices[0] = {1, 2, 3};
        mi_complex.vertices[1] = {2, 0, 3};
        mi_complex.vertices[2] = {0, 1, 3};

        mi_complex.edges.resize(3);
        mi_complex.edges[0].vertices = {1, 2};
        mi_complex.edges[1].vertices = {2, 0};
        mi_complex.edges[2].vertices = {0, 1};

        mi_complex.edges[0].positive_material_label = 0;
        mi_complex.edges[0].negative_material_label = 3;
        mi_complex.edges[1].positive_material_label = 1;
        mi_complex.edges[1].negative_material_label = 3;
        mi_complex.edges[2].positive_material_label = 2;
        mi_complex.edges[2].negative_material_label = 3;

        mi_complex.faces.resize(1);
        mi_complex.faces[0].material_label = 3;
        mi_complex.faces[0].edges = {0, 1, 2};
    } else {
        mi_complex.vertices[0] = {1, 2, 3, 4};
        mi_complex.vertices[1] = {0, 2, 3, 4};
        mi_complex.vertices[2] = {0, 1, 3, 4};
        mi_complex.vertices[3] = {0, 1, 2, 4};

        mi_complex.edges.resize(6);
        mi_complex.edges[0].vertices = {0, 1};
        mi_complex.edges[1].vertices = {0, 2};
        mi_complex.edges[2].vertices = {0, 3};
        mi_complex.edges[3].vertices = {1, 2};
        mi_complex.edges[4].vertices = {1, 3};
        mi_complex.edges[5].vertices = {2, 3};

        mi_complex.edges[0].supporting_materials = {2, 3, 4};
        mi_complex.edges[1].supporting_materials = {1, 3, 4};
        mi_complex.edges[2].supporting_materials = {1, 2, 4};
        mi_complex.edges[3].supporting_materials = {0, 3, 4};
        mi_complex.edges[4].supporting_materials = {0, 2, 4};
        mi_complex.edges[5].supporting_materials = {0, 1, 4};

        mi_complex.faces.resize(4);
        mi_complex.faces[0].edges = {3, 5, 4};
        mi_complex.faces[1].edges = {1, 2, 5};
        mi_complex.faces[2].edges = {0, 4, 2};
        mi_complex.faces[3].edges = {0, 1, 3};

        mi_complex.faces[0].positive_material_label = 0;
        mi_complex.faces[0].negative_material_label = 4;
        mi_complex.faces[1].positive_material_label = 1;
        mi_complex.faces[1].negative_material_label = 4;
        mi_complex.faces[2].positive_material_label = 2;
        mi_complex.faces[2].negative_material_label = 4;
        mi_complex.faces[3].positive_material_label = 3;
        mi_complex.faces[3].negative_material_label = 4;

        mi_complex.cells.resize(1);
        mi_complex.cells[0].material_label = 4;
        mi_complex.cells[0].faces = {0, 1, 2, 3};
    }
    return mi_complex;
}

} // namespace simplicial_arrangement
