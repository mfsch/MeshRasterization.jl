# Compute Rasterized Representations of Polygon Meshes

MeshRasterization.jl provides functionality to compute a number of geometric quantities of
a mesh for a large number of points, usually arranged in a Cartesian grid with regular or
irregular grid spacing.

The following (related) quantities can currently be computed:

- closest point on mesh boundary
- unsigned and signed distance to the closest point
- direction of the closest point

The package implements a version of the Characteristic/Scan Conversion (CSC) algorithm that
was originally introduced by Mauch (2000). The implementation mostly follows the description
of the algorithm by Roosing et al. (2019), which relies on angle-weighted pseudonormals
(Bærentzen & Aanæs, 2005) for more robust results.

The package also implements a version of the scan conversion algorithm originally described
by Pineda (1988) to efficiently find grid points that lie within a convex polyhedron.

The implementation relies on the [Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl)
for geometric types and some of the algorithms.


## Rasterization

The functions `distance`, `signed_distance`, `closest_point`, and `direction` compute the
respective quantities for a given mesh and a list of raster points. Please refer to the
[docstrings](https://docs.julialang.org/en/v1/manual/documentation/) of those functions for
details.

The function `rasterize` can be used to compute multiple quantities in a single pass, and
the function `rasterize!` writes the results to pre-allocated arrays. The latter also allows for more fine-grained control over the return types.


## Scan Conversion

The function `scan(geometry, points[, method])` from the submodule `ScanConversion` produces
an iterator over the indices of `points` that lie within a Meshes.jl `geometry`. For convex
polyhedra, the function can use the `EdgeFunctionScan` method, which is much more efficient
and can handle points on the boundary of the geometry in a robust way. The function does not
attempt to detect convex polyhedra except for a few trivial geometries, so the
`EdgeFunctionScan` method usually has to be specified explicitly.

The `points` can be a `CartesianGrid` from Meshes.jl, a tuple of coordinate axes, or
a collection of `Point`s from Meshes.jl. For grid inputs, the location of cell centroids is
used. The iterator returns the `CartesianIndex` of each point inside the `geometry`, or
elements of `eachindex(points)` for collection of `Point`s.


## Example

Computing a signed distance field:

```julia
using Meshes, MeshRasterization
points = [(1,1,1),(2,1,1),(1,2,1),(1,1,2)]
connec = [(1,2,3),(1,3,4),(1,4,2),(4,3,2)]
mesh = SimpleMesh(points, connect.(connec))
grid = CartesianGrid((0.,0.,0.), (3.,3.,3.), dims=(10,10,10))
sdf = signed_distance(mesh, grid, dmax=0.5)
```

Computing several fields at once:

```
rasterize((:signed_distance, :direction, :closest_point), mesh, grid, dmax=1)
```

Finding indices inside a geometry:

```julia
julia> using MeshRasterization.ScanConversion

julia> geom = Ball((1, 1, 1.5), 1)
Ball{3,Float64}(Point(1.5, 1.0, 1.0), 1.0))

julia> grid = CartesianGrid(1, 2, 3)
1×2×3 CartesianGrid{3,Float64}
  minimum: Point(0.0, 0.0, 0.0)
  maximum: Point(1.0, 2.0, 3.0)
  spacing: (1.0, 1.0, 1.0)

julia> axes = (0.5:1:1, 0.5:1:2, 0.5:1:3)
(0.5:1.0:0.5, 0.5:1.0:1.5, 0.5:1.0:2.5)

julia> collect(scan(geom, grid))
2-element Vector{CartesianIndex{3}}:
 CartesianIndex(1, 1, 2)
 CartesianIndex(1, 2, 2)

julia> collect(scan(geom, axes))
2-element Vector{CartesianIndex{3}}:
 CartesianIndex(1, 1, 2)
 CartesianIndex(1, 2, 2)
```


## References

- Bærentzen J. A. and Aanæs H. (2005). Signed distance computation using the angle weighted
  pseudonormal. *IEEE Transactions on Visualization and Computer Graphics* 11 (3), 243–253.
  [doi:10.1109/TVCG.2005.49](https://doi.org/10.1109/TVCG.2005.49).
- Mauch, S. (2000). A fast algorithm for computing the closest point and distance transform.
  *Caltech ASCI Technical Report* 077.
- Pineda, J. (1988). A parallel algorithm for polygon rasterization. SIGGRAPH Computer Graphics 22 (4), 17–20. [doi:10.1145/378456.378457](https://doi.org/10.1145/378456.378457).
- Roosing, A., Strickson, O. and Nikiforakis, N. (2019). Fast distance fields for fluid
  dynamics mesh generation on graphics hardware. *Communications in Computational Physics*
  26 (3), 654-680. [doi:10.4208/cicp.OA-2018-013](https://doi.org/10.4208/cicp.OA-2018-013).
