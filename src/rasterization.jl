module Rasterization

import LinearAlgebra: dot, norm
import Meshes: Point, Vec, normal, CartesianGrid, centroid,
    SimpleMesh, vertices, topology, nelements, element,
    Coboundary, indices, materialize

# ----------------
# Shared Interface
# ----------------

abstract type RasterizationMethod end
function rasterize! end
function rasterize end


# ----------------------------
# Raster Point Representations
# ----------------------------

pointdims(points::Tuple) = length.(points)
pointdims(points) = size(points)

getpoint(points, ind) = points[ind]
getpoint(points::Tuple, ind) = Point(getindex.(points, Tuple(ind)))
getpoint(grid::CartesianGrid, ind) = centroid(grid[ind])

pointindices(points) = eachindex(points)
pointindices(points::Tuple) = CartesianIndices(eachindex.(points))
pointindices(grid::CartesianGrid) = CartesianIndices(size(grid))


# -----------------
# Output Management
# -----------------

infinity(::Type{Float64}) = Inf
infinity(::Type{Float32}) = Inf32

function current_distance(data, ind)
    if haskey(data, :distance)
        data.distance[ind]
    elseif haskey(data, :signed_distance)
        abs(data.signed_distance)[ind]
    else
        error("Data does not contain any field for distance")
    end
end

function reset_distance!(data)
    for key in (:distance, :signed_distance)
        haskey(data, key) || continue
        values = getproperty(data, key)
        fill!(values, infinity(eltype(values)))
    end
end


# ----------------------
# Weighted Pseudonormals
# ----------------------

function vertex_orientation(mesh::SimpleMesh{3}, ind)
    sum(Coboundary{0,2}(topology(mesh))(ind)) do el
        conn = element(topology(mesh), el)
        vert = vertices(mesh)

        # compute normal vector of triangle
        n = normal(materialize(conn, vert))

        # compute angle at vertex
        v1, v2, v3 = indices(conn)
        @assert ind in (v1, v2, v3)
        n1 = v1 == ind ? v2 : v1 # first neighbor
        n2 = v3 == ind ? v2 : v3 # second neighbor
        vec1 = vert[n1] - vert[ind]
        vec2 = vert[n2] - vert[ind]
        α = acos(dot(vec1, vec2) / (norm(vec1) * norm(vec2)))
        α * n # result is not normalized!
    end
end

function vertex_orientation(mesh::SimpleMesh{2,T}, ind) where {T}
    elem4vert = prepare_elem4vert(topology(mesh))
    dir = zero(Vec{2,T})
    for el in elem4vert[ind]
        isnothing(el) && continue
        dir += normal(element(mesh, el))
    end
    dir
end

function prepare_elem4vert(topo)
    elem4vert = Dict{Int,NTuple{2,Union{Int,Nothing}}}()
    for el in 1:nelements(topo)
        v1, v2 = indices(element(topo, el))

        # update entries for vertex before current edge
        before1, after1 = get(elem4vert, v1, (nothing, nothing))
        isnothing(after1) || error("Vertex $v1 is followed by two edges")
        elem4vert[v1] = (before1, el)

        # update entries for vertex after current edge
        before2, after2 = get(elem4vert, v2, (nothing, nothing))
        isnothing(before2) || error("Vertex $v2 is preceeded by two edges")
        elem4vert[v2] = (el, after2)
    end
    elem4vert
end

end
