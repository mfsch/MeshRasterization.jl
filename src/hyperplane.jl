module Hyperplanes

using LinearAlgebra: cross
using Meshes: Point, Vec, Ngon, Polyhedron, Polytope, Primitive,
              vertices, facets

import Meshes: normal

using ..Prisms: Prism

struct Hyperplane{Dim,T} <: Primitive{Dim,T}
    p::Point{Dim,T}
    n::Vec{Dim,T}
end

normal(hp::Hyperplane) = hp.n

"""
hyperplanes(polytope)

Return the hyperplanes in which the boundaries of a polytope live, oriented such that the plane normals point out of the polytope.
"""
function hyperplanes(pt::Polytope{Dim,Dim,T}) where {Dim,T}
    error("Not implemented")
end

function hyperplanes(pg::Ngon{N,2,T}) where {N,T}
    v = vertices(pg)
    ntuple(N) do i
        p = v[i]
        dp = p - (i == 1 ? v[end] : v[i-1])
        n = Vec(dp[2], -dp[1])
        Hyperplane(p, n)
    end
end

function hyperplanes(ph::Polyhedron{3,T}) where {T}
    # TODO: this only works correctly for some polyhedra
    Iterators.map(facets(ph)) do facet
        xs = vertices(facet)
        p::Point{3,T} = xs[1]
        dp1::Vec{3,T} = xs[2] - p
        dp2::Vec{3,T} = xs[3] - p
        Hyperplane(p, cross(dp1, dp2))
    end
end

function hyperplanes(pr::Prism{3,T}) where {T}
    v1 = pr.base[2] - pr.base[1]
    v2 = pr.base[3] - pr.base[2]
    v3 = pr.base[1] - pr.base[3]
    v4 = pr.side
    vs = (v1, v2, v3, v4)
    ntuple(5) do i
        p = isodd(i) ? pr.base[1] : pr.base[2] + pr.side
        dp1 = vs[(1,2,3,1,1)[i]]
        dp2 = vs[(4,4,4,2,3)[i]]
        Hyperplane(p, cross(dp1, dp2))
    end
end

end
