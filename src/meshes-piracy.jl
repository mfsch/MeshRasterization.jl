"""
    MeshesPiracy

New method implementations for functions defined either in Meshes.jl or elsewhere that act
on types of Meshes.jl. Since this is a case of [“type
piracy”](https://docs.julialang.org/en/v1/manual/style-guide/#Avoid-type-piracy), those
methods should ideally be moved to Meshes.jl.
"""
module MeshesPiracy

using LinearAlgebra: norm
using Meshes: Point, Vec, Segment, Ngon, Polytope, SimpleTopology,
              HalfEdgeTopology, CartesianGrid, vertices, indices, Connectivity,
              connect, coordinates

import Base: zero, convert, contains
import Meshes: nfacets, facets, nfaces, faces, element, normal, centroid

# allow cartesian indices for centroid function
function centroid(g::CartesianGrid{Dim}, ind::CartesianIndex) where {Dim}
  neworigin = coordinates(g.origin) .+ g.spacing ./ 2
  Point(ntuple(i -> neworigin[i] + (ind[i] - g.offset[i])*g.spacing[i], Dim))
end

# compute normals for other dimensions of polygons etc.
function normal(ngon::Ngon{2,2})
    a, b = vertices(ngon)
    v = b - a
    n = Vec(v[2], -v[1])
    n / norm(n)
end
normal(s::Segment{2}) = normal(Ngon(vertices(s)))

# define the zero element for `Point` as the origin
zero(::T) where {T<:Point} = zero(T)
zero(::Type{Point{Dim,T}}) where {Dim,T} = Point((zero(T) for _ in 1:Dim)...)

# avoid duplicating work of half-edge topology initialization
convert(::Type{HalfEdgeTopology}, topo::HalfEdgeTopology) = topo

# allow using CartesianIndex for CartesianGrid
function element(g::CartesianGrid{Dim}, ind::CartesianIndex{Dim}) where {Dim}
    element(g, LinearIndices(size(g))[ind])
end

# define the faces not only for the boundary of a polytope but also for the polytope itself
# to avoid constructing the boundary where this is not needed
function faces(p::Polytope, rank)
    if rank == 0
        vertices(p)
    elseif rank == 1
        edges(p)
    elseif rank == paramdim(p) - 1
        facets(p)
    else
        faces(boundary(p), 1)
    end
end

function faces(c::Connectivity{<:Ngon{N},N}, rank) where N
    if rank == 0
        indices(c)
    elseif rank == 1
        ntuple(N) do i
            connect((i, i == N ? 1 : i+1))
        end
    elseif rank == 2
        conn
    else
        error("Rank '$rank' cannot be larger than parametric dimension '2'")
    end
end

function nfaces(p::Polytope, rank)
    if rank == 0
        nvertices(p)
    elseif rank == 1
        nedges(p)
    elseif rank == paramdim(p) - 1
        nfacets(p)
    else
        nfaces(boundary(p), rank)
    end
end

nfacets(p::Polytope) = nfaces(boundary(p), paramdim(p)-1)

facets(p::Polytope) = faces(boundary(p), paramdim(p)-1)


end
