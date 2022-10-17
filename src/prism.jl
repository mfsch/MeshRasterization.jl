module Prisms

import LinearAlgebra: dot, cross

using Meshes: Point, Vec, Polyhedron, SimpleMesh, connect, materialize
import Meshes: nvertices, vertices, nfacets, facets, boundary

struct Prism{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Polyhedron{Dim,T}
    base::V # should be oriented such that the normal and the side are colinear
    side::Vec{Dim,T}

    function Prism(base::V, side, fix=true) where {Dim,T,V<:AbstractVector{Point{Dim,T}}}
        if fix && sign(dot(cross(base[2]-base[1], base[3]-base[1]), side)) != 1
            base = reverse(base)
        end
        new{Dim,T,V}(base, side)
    end
end


function connections(p::Prism, rank)
    if rank == 0
        ((1,), (2,), (3,), (4,), (5,), (6,))
    elseif rank == 1
        ((1,2), (2,3), (3,1), (1,4), (2,5), (3,6), (4,5), (5,6), (6,4))
    elseif rank == 2
        ((3,2,1),(5,4,1,2),(3,6,5,2),(1,4,6,3),(5,6,4))
    elseif rank == 3
        ((1, 2, 3, 4, 5, 6))
    else
        error("Invalid rank '$rank' for prism")
    end
end

# NOTE: creating a boundary allocates a vector for the connections
boundary(p::Prism) = SimpleMesh(vertices(p), collect(connect.(connections(p, 2))))

nvertices(p::Prism) = 6
function vertices(p::Prism)
    Vec(p.base[1], p.base[2], p.base[3], p.base[1] + p.side, p.base[2] + p.side,
        p.base[3] + p.side)
end

nfacets(p::Prism) = 5

function facets(p::Prism)
    v = vertices(p)
    c = connections(p, 2)
    (materialize(connect(c), v) for c in c)
end

end
