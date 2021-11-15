# Simplicial Arrangement

## Overview

A simplicial arrangement is defined as the arrangement over a 2 or 3 dimensional
simplex cut by a set of hyperplanes.

## Build

```sh
mkdir build
cd build

# Option 1: just build the library.
cmake ..
make

# Option 2: build both library and unit tests.
cmake .. -DSIMPLICIAL_ARRANGEMENT_UNIT_TESTS=On
make
./simplicial_arrangement_tests           # to run unit tests.
./simplicial_arrangement_tests benchmark # to run benchmark.
```

## Quick start
This library provides a single function `compute_arrangement` and the
corresponding arrangement data structures.

```c++
#include <simplicial_arrangement/simplicial_arrangement.h>

using namespace simplicial_arrangement;
using Scalar = double;

std::vector<Plane<Scalar, 3>> cut_planes;
// ... populate cut_planes ...

auto arrangement = compute_arrangement(cut_planes);
```
The output `arrangement` object represents the arrangement induced by the cut
planes within a simplex.  It contains information of vertices, faces and cells
as well as the cell adjacency and coplanar planes. Please see [data
structures](#data-structures) section for more details.

## Data structures

### Plane
A plane is represented implicitly as the coefficients of the barycentric plane
equation:

```c++
using namespace simpicial_arrangement;

// For 2D, the plane (f0, f1, f2) means
// f0 * b0 + f1 * b1 + f2 * b2 = 0
// where b0, b1, b2 are variables representing barycentric coordinates.
Plane<double, 2> plane_2d{f0, f1, f2};

// For 3D, the plane (f0, f1, f2, f3) means
// f0 * b0 + f1 * b1 + f2 * b2 + f3 * b3 = 0
// where b0, b1, b2, b3 are variables representing barycentric coordinates.
Plane<double, 3> plane_3d{f0, f1, f2, f3};
```

### Point
A point is represented implicitly as the intersection of `dim` planes, where
`dim` is either 2 or 3:

```c++
Point<2> point_2d{0, 1};    ///< 2D intersection of plane 0 and 1.
Point<3> point_3d{0, 1, 2}; ///< 3D intersection of plane 0, 1, and 2.
```

### Arrangement
The `Arrangement<DIM>` structure contains information about the vertices, faces,
cells and unique planes.  Together, they provide a flattened representation of
the arrangement induced by the input planes within a `DIM`-dimensional simplex.

Note that within `Arrangement<DIM>` structure, we use the convention that the
first `DIM+1` planes are the boundary planes of the simplex.  The `i`th user
provided cut plane will have index `DIM+1+i`.

#### Vertices
The `vertices` of an arrangement include all intersection points as well as the
corners of the simplex.  They are represented as an array of `Point<DIM>`s.

**Guarantee**: the arrangement vertices contain no duplicates.

#### Faces
The `faces` of an arrangement represents a set of (`DIM-1`)-dimensional
polytopes induced by the cut planes within the simplex.  Each face is
represented by the inner `Face` data structure.

In 2D, a face is simply an edge.
In 3D, a face is a convex polygon.

* `Face::vertices` represents the end points of the 2D edge or the boundary
  vertex loop of a 3D polygon.
* `Face::supporting_plane` is a plane that contains the face
  <sup>[1](#coplanar_planes)</sup>.
* `Face::positive_cell` and `Face::negative_cell` are cell indices to the
  corresponding cell on the positive and negative side of the supporting plane
  <sup>[2](#boundary_face)</sup>.

**Guarantees**:
1. In 2D, `Face::vertices` are ordered such that the positive side
of the `Face::supporting_plane` is on the right side of the edge.
2. In 3D, `Face::vertices` are in counterclockwise order when viewed from the positive
side of the supporting plane.
3. There is no duplicate faces regardless of the orientation.  I.e. no 2 faces
   are coplanar.

<a name="coplanar_planes">**Note<sup>1</sup>**</a>:
If there are multiple planes containing a face, `Face::supporting_plane` can be
any of these planes.  Please see the [unique planes section](#unique-planes) for
data that allow one to recover all planes that are coplanar with this face.

<a name="boundary_face">**Note<sup>2</sup>**</a>:
If a face is on the boundary of the simplex, the value of
`Face::positive_cell` or `Face::negative_cell` may be `None`.

#### Cells
The `cells` of an arrangement represents a set of `DIM`-dimensional convex
polytopes.  Each cell is represented by the inner `Cell` data structure.

* `Cell::faces` is a list of faces that form the boundary of this cell.  In 2D,
  they in counterclockwise order.  In 3D, they are unordered.
* `Cell::face_orientations` is a list of orientations of the faces relative to
  the cell.  I.e. `face_orientations[i] == true` means the cell is on the
  positive side of the `i`th face.
* `Cell::plane_orientations` is a list of orientations of the cell with respect
  to all input planes.  I.e. `plane_orientation[i] == true` means the cell is on
  the positive side of the `i`th plane.

#### Unique planes
If the input cut planes contain coplanar planes, their information is stored in
three lists:

```c++
// To find the unique plane index of plane i
size_t uid = arrangement.unique_plane_indices[i];

// To get the set of coplanar planes corresponding to uid
auto coplanar_planes = arrangement.unique_planes[uid];

// Assume plane j and k are coplanar, to check the relative orientation of
// jth and kth planes:
if (arrangement.plane_orientations[j] == arrangement.plane_orientations[k]) {
    // plane j and k have the same orientation.
}
```

