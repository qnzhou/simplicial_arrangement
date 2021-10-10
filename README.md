# Simplicial Arrangement

## Overview

A simplicial arrangement is defined as the arrangement over a 2 or 3 dimensional
simplex cut by a set of hyperplanes.

In this work, a hyperplane is represented implicitly by the coefficients of the
Barycentric plane equation.  A point is represented indirectly as the
intersection of `dim` hyperplanes, where `dim` is either 2 or 3.

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
./simplicial_arrangement_tests
```

## Usage

```c++
#include <simplicial_arrangement/simplicial_arrangement.h>

// 3D with double as scalar.
std::vector<simplicial_arrangement::Plane<double, 3>> planes;
// Populate planes.
auto arrangement = simplicial_arrangement::compute_arrangement(planes);
```

The output `arrangement` object contains information of vertices, faces and
cells as well as the cell adjacency information.  Please see
[simplicial_arrangement.h](include/simplicial_arrangement/simplicial_arrangement.h)
for more details.

## Data structures

### Plane
A plane is represented implicitly as the coefficents of the barycentric plane
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
Point<2> p_2d{0, 1};    ///< 2D intersection of plane 0 and 1.
Point<3> p_3d{0, 1, 2}; ///< 3D intersection of plane 0, 1, and 2.
```

### Arrangement
The `Arrangement<DIM>` structure contains information about the vertices, faces,
cells and unique planes.  Together, they provide a flattened representation of
the arrangement induced by the input planes within a `DIM`-dimensional simplex.

Note that within `Arrangement<DIM>` structure, we use the convention that the
first `DIM+1` planes are the boundary planes of the simplex.  The `i`th user
provided cut plane will be off index `DIM+1+i`.

#### Vertices
The `vertices` of the arrangement include all intersection points as well as the
corners of the simplex.  They are represented as an array of `Point<DIM>`s.

**Guarantee**: the arrangement vertices contain no duplicates.

#### Faces
The `faces` of the arrangement represents a set of (`DIM-1`)-dimensional
polytopes induced by the cut planes within the simplex.  Each face is
represented by the inner `Face` data structure.

##### 2D
In 2D, a face is simply an edge.

* `Face::vertices` represents the end points of the 2D edge.
* `Face::supporting_plane` is a line that contains the 2D edge.
* `Face::positive_cell` and `Face::negative_cell` are cell indices to the
  corresponding cell on the positive and negative side of the supporting plane
  <sup>[1](#boundary_face)</sup>.

**Guarantee**: `Face::vertices` are ordered such that the positive side of the
`Face::supporting_plane` is on the right side of the edge.

##### 3D

In 3D, a face is a convex polygon.

* `Face::vertices` represents the vertex loop that forms the polygon boundary.
* `Face::supporting_plane` is a plane that contains the polygon.
* `Face::positive_cell` and `Face::negative_cell` are cell indices to the
  corresponding cell on the positive and negative side of the supporting plane
  <sup>[1](#boundary_face)</sup>.

**Guarantee**: `Face::vertices` are in counterclockwise order when viewed from
the positive side of the supporting plane.

<a name="boundary_face">**Note<sup>1</sup>**</a>:
If a face is on the boundary of the simplex, the value of
`Face::positive_cell` or `Face::negative_cell` may be `None`.

#### Cells

#### Unique planes

