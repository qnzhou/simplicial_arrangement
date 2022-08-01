#Simplicial Arrangement

## Overview

This repo contains the simplicial kernel implementations for our Siggraph 2022
paper "[Robust Computation of Implicit Surface Networks for Piecewise Linear
Functions]".  Specifically, it robustly computes the arrangement or material
interface induced by a set of implicit functions within a single simplex in 2D
and 3D.

## Build

```sh
mkdir build
cd build

#Option 1 : just build the library.
cmake ..
make

#Option 2 : build both library and unit tests.
cmake .. -DSIMPLICIAL_ARRANGEMENT_UNIT_TESTS=On
make
./simplicial_arrangement_tests           # to run unit tests.
./simplicial_arrangement_tests benchmark # to run benchmark.
```

## Quick start

This library provides a minimalistic single function interface,
`compute_arrangement` and `compute_material_interface`, that is self-explanatory.

```c++
#include <simplicial_arrangement/simplicial_arrangement.h>

using namespace simplicial_arrangement;
using Scalar = double;

// For arrangement computation.
std::vector<Plane<Scalar, 3>> cut_planes;
// ... populate cut_planes ...
auto arrangement = compute_arrangement(cut_planes);

// For material interface computation.
std::vector<Material<Scalar, 3>> materials;
// ... populate materials ...
auto material_interface = compute_material_interface(materials);
```

The output `arrangement` and `material_interface` objects represent the space
partition induced by arrangement and materials.  Both data structures contain
information of vertices, faces and cells as well as the cell adjacency. Please
see data structure sections ([arrangement](#Arrangement-data-structures),
[material interface](#material-interface-data-structures)) for more details.

## Arrangement data structures

### Plane

The zero crossing of an implicit function within a linear simplicial element is
a plane. A plane is represented implicitly as the coefficients of the
barycentric plane equation:

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
`dim` is either 2 or 3 :

```c++
Point<2> point_2d{0, 1}; ///< 2D intersection of plane 0 and 1.
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
if (arrangement.unique_plane_orientations[j] == arrangement.unique_plane_orientations[k]) {
    // plane j and k have the same orientation.
}
```

## Material interface data structures

### Material

A material is represented as an implicit function's evaluation at the vertices
of a simplex.

```c++
using namespace simplicial_arrangement;

// For 2D material, let f0, f1, f2 be the implicit function values at vertices of a triangle.
Material<double, 2> m_2d{f0, f1, f2};

// For 3D material, let f0, f1, f2, f3 be the implicit function values at vertices of a tet.
Material<double, 2> m_3d{f0, f1, f2};
```

### Joint

A joint is a point where at least `dim+1` materials evaluates to the same value
(`dim` is either 2 or 3).

```c++
Joint<2> joint_2d{0, 1, 2}; ///< The 2D point where implicit functions 0, 1, 2 are the same.
Joint<3> joint_3d{0, 1, 2, 3}; ///< The 3D point where implicit functions 0, 1, 2, 3 are the same.
```

### MaterialInterface

The `MaterialInterface<DIM>` structure contains information about the vertices,
faces, cells and either adjacency.  Together, they define a space partition
induced by a set of materials within a `DIM`-dimensional simplex.

Note that we create pseudo material to denote simplex faces.  The first `DIM+1`
materials corresponds to these pseudo materials.  The `i`th user provided
material will have index `DIM+i+1` in the output.

#### Vertices
The `vertices` of an arrangement include all material joints as well as corners
of the simplex.  They are represented as an array of `Joint<DIM>`.

**Guarantee**: the material interface vertices contains no duplicates.

#### Faces
The `faces` of a material interface is a set of (`DIM=1`)-dimensional polytopes
where 2 of the materials evaluates to the same value.  Each face contains the
following fields:

In 2D, a face is simply an edge.
In 3D, a face is a convex polygon.

* `Face::vertices` presents the end points of a 2D edge or the boundary loop of a
  3D polygon.
* `Face::positive_material_label` is the material id on the positive side of the
  face.
* `Face::negative_material_label` is the material id on the negative side of the
  face.

**Guarantees**:
1. In 2D, `Face::vertices` are ordered such that the positive side
of the `Face::supporting_plane` is on the right side of the edge.
2. In 3D, `Face::vertices` are in counterclockwise order when viewed from the positive
material side.
3. There is no duplicate faces regardless of the orientation.  I.e. no 2 faces
   are coplanar.

#### Cells
The `cells` of a material interface represents a `DIM`-dimensional convex
polytope where a single material dominates (i.e. its values is larger than all
other materials for any point within the cell).  Each cell is represented by the
inner `Cell` data structure.

* `Cell::faces` is a list of faces that form the boundary of the cell.  In 2D,
  they are in counterclockwise order.  In 3D, they are unordered.
* `Cell::material_label` is the id of the dominant material in this cell.

#### Unique materials

If the input contains materials that have the same evaluation on the simplex
vertices, they are considered as duplicate materials.  Their information is
stored in the following way:

* `unique_material_indices` is a mapping from the input material id to the
  unique material id.  It is empty if all materials are unique.
* `unique_materials` is a list of unique material groups.  I.e. all material
  listed in `unique_materials[i]` are the same and they share the same unique
  material id `i`.  This list is empty if all materials are unique.


[Robust Computation of Implicit Surface Networks for Piecewise Linear Functions]: https://duxingyi-charles.github.io/publication/robust-computation-of-implicit-surface-networks-for-piecewise-linear-functions/
