#pragma once

#include "ARComplex.h"

namespace simplicial_arrangement {

template <int DIM>
ARComplex<DIM> initialize_simplicial_ar_complex(size_t num_planes)
{
    ARComplex<DIM> ar_complex;
    ar_complex.vertices.resize(DIM + 1);

    if constexpr (DIM == 2) {
        ar_complex.vertices[0] = {1, 2};
        ar_complex.vertices[1] = {2, 0};
        ar_complex.vertices[2] = {0, 1};

        ar_complex.edges.resize(3);
        ar_complex.edges[0].vertices = {2, 1};
        ar_complex.edges[1].vertices = {0, 2};
        ar_complex.edges[2].vertices = {1, 0};
        ar_complex.edges[0].supporting_plane = 0;
        ar_complex.edges[1].supporting_plane = 1;
        ar_complex.edges[2].supporting_plane = 2;
        ar_complex.edges[0].positive_face = 0;
        ar_complex.edges[0].negative_face = INVALID;
        ar_complex.edges[1].positive_face = 0;
        ar_complex.edges[1].negative_face = INVALID;
        ar_complex.edges[2].positive_face = 0;
        ar_complex.edges[2].negative_face = INVALID;

        ar_complex.faces.resize(1);
        ar_complex.faces[0].edges = {0, 1, 2};
        ar_complex.faces[0].signs = {true, true, true};
        ar_complex.faces[0].signs.resize(num_planes, false);
    } else {
        ar_complex.vertices[0] = {1, 2, 3};
        ar_complex.vertices[1] = {2, 3, 0};
        ar_complex.vertices[2] = {3, 0, 1};
        ar_complex.vertices[3] = {0, 1, 2};

        ar_complex.edges.resize(6);
        ar_complex.edges[0].vertices = {0, 1};
        ar_complex.edges[1].vertices = {0, 2};
        ar_complex.edges[2].vertices = {0, 3};
        ar_complex.edges[3].vertices = {1, 2};
        ar_complex.edges[4].vertices = {1, 3};
        ar_complex.edges[5].vertices = {2, 3};
        ar_complex.edges[0].supporting_planes = {2, 3};
        ar_complex.edges[1].supporting_planes = {1, 3};
        ar_complex.edges[2].supporting_planes = {1, 2};
        ar_complex.edges[3].supporting_planes = {0, 3};
        ar_complex.edges[4].supporting_planes = {0, 2};
        ar_complex.edges[5].supporting_planes = {0, 1};

        ar_complex.faces.resize(4);
        ar_complex.faces[0].edges = {5, 3, 4};
        ar_complex.faces[1].edges = {2, 1, 5};
        ar_complex.faces[2].edges = {4, 0, 2};
        ar_complex.faces[3].edges = {1, 0, 3};
        ar_complex.faces[0].supporting_plane = 0;
        ar_complex.faces[1].supporting_plane = 1;
        ar_complex.faces[2].supporting_plane = 2;
        ar_complex.faces[3].supporting_plane = 3;
        ar_complex.faces[0].positive_cell = 0;
        ar_complex.faces[0].negative_cell = INVALID;
        ar_complex.faces[1].positive_cell = 0;
        ar_complex.faces[1].negative_cell = INVALID;
        ar_complex.faces[2].positive_cell = 0;
        ar_complex.faces[2].negative_cell = INVALID;
        ar_complex.faces[3].positive_cell = 0;
        ar_complex.faces[3].negative_cell = INVALID;

        ar_complex.cells.resize(1);
        ar_complex.cells[0].faces = {0, 1, 2, 3};
        ar_complex.cells[0].signs = {true, true, true, true};
        ar_complex.cells[0].signs.resize(num_planes, false);
    }
    return ar_complex;
}

} // namespace simplicial_arrangement
