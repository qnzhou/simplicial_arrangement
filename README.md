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
