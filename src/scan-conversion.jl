module ScanConversion

using LinearAlgebra: dot
using Meshes: Geometry, CartesianGrid, Box, Vec, Polytope,
              boundingbox, spacing

using ..Hyperplanes: hyperplanes

export SimpleScan, EdgeFunctionScan, scan

"""
    ScanMethod

A method for finding the subset of regularly spaced points that are inside a geometry.
"""
abstract type ScanMethod end

"""
    SimpleScan

Use the `in` function to determine which points are inside a geometry.
"""
struct SimpleScan <: ScanMethod end

"""
    EdgeFunctionScan

Use the concept of edge functions as described by Pineda (1988) to efficiently determine
which points are inside a geometry.

!!! warn
    This approach only works for convex polhyedra.

## References

* [Juan Pineda. 1988. A parallel algorithm for polygon rasterization. SIGGRAPH Comput.
  Graph. 22, 4 (Aug. 1988), 17–20.](https://doi.org/10.1145/378456.378457)
"""
struct EdgeFunctionScan{T} <: ScanMethod
    tol::T
    EdgeFunctionScan(tol = 1e-9) = new{typeof(tol)}(tol)
end

# use SimpleScan by default
scanmethod(geom, raster) = SimpleScan()


"""
    scan(geometry, raster[, method])

Find the subset of the points described by the raster that are inside the geometry.
"""
function scan end

scan(geom, raster) = scan(geom, raster, scanmethod(geom, raster))



function scan(box::Box{Dim,T},
        grid::CartesianGrid{Dim,T}) where {Dim,T}
    dx = spacing(grid)
    imin = (minimum(box) - minimum(grid)) ./ dx .+ 1 |>
    i -> floor.(Int, i) |>
    i -> max.(i, 1)
    imax = (maximum(box) - minimum(grid)) ./ dx |>
    i -> ceil.(Int, i) |>
    i -> min.(i, size(grid))
    inds = ntuple(i -> imin[i]:imax[i], Dim)
    CartesianIndices(inds)
end

# test if point in polygon for all points in bounding box
function scan(geom::Geometry{Dim,T},
        grid::CartesianGrid{Dim,T},
        ::SimpleScan) where {Dim,T}

    Δgd = spacing(grid) |> Vec
    xref = minimum(grid) - 0.5 * Δgd

    inds = scan(boundingbox(geom), grid)
    Iterators.filter(inds) do i
        pt = xref + Δgd .* i.I
        pt ∈ geom
    end
end

function scan(polytope::Polytope{Dim,Dim,T},
        grid::CartesianGrid{Dim,T},
       method::EdgeFunctionScan) where {Dim,T}
    E = edgefunctions(polytope, grid)
    tol = convert(T, method.tol)
    inds = scan(boundingbox(polytope), grid)
    Iterators.filter(inds) do ind
        all(apply_edgefunction(Ei, ind) < tol for Ei in E)
    end
end

function edgefunctions(polytope::Polytope{Dim,Dim,T},
        grid::CartesianGrid{Dim,T},
        pos = T(0.5)) where {Dim,T}
    dx = spacing(grid) |> Vec
    xref = minimum(grid) + (pos .- 1) .* dx
    map(hyperplanes(polytope)) do hp
        # as the facets should have outward-pointing normals, the edge functions
        # as defined here should be negative for points inside the polyhedron
        dE = hp.n .* dx
        Eref = dot(xref - hp.p, hp.n)
        Vec(dE..., Eref)
    end
end

# explicit multiplication appears to run faster than dot product
@inline apply_edgefunction(Ei::Vec{3}, ind::CartesianIndex{2}) =
Ei[1] * ind[1] + Ei[2] * ind[2] + Ei[3]
@inline apply_edgefunction(Ei::Vec{4}, ind::CartesianIndex{3}) =
Ei[1] * ind[1] + Ei[2] * ind[2] + Ei[3] * ind[3] + Ei[4]

end
