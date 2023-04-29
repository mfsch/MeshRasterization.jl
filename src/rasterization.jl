module Rasterization

using LinearAlgebra: dot, norm
using Meshes: Vec, Point, SimpleMesh, HalfEdgeTopology, normal, vertices,
              topology, nelements, element, Coboundary, indices, materialize
using ..Rasters: Raster, dimensions

export closest_point, direction, distance, signed_distance, rasterize, rasterize!

abstract type RasterizationMethod end
function rasterize! end
default_method() = nothing

function closest_point(mesh, points)
    rasterize((:closest_point,), mesh, points).closest_point
end

function direction(mesh, points, args...)
    rasterize((:direction,), mesh, points, args...).direction
end

function distance(mesh, points, args...)
    rasterize((:distance,), mesh, points, args...).distance
end

function signed_distance(mesh, points, args...)
    rasterize((:signed_distance,), mesh, points, args...).signed_distance
end

"""
    rasterize(fields, mesh, points)

Compute rasterized representations of a `mesh` at a collection of `points`.

# Arguments

- `fields`: a list of symbols that define the raster data included in the output
  - `:closest_point`: the closest point on the mesh boundary (`Point`)
  - `:direction`: unit vector in the direction of the closest point (`Vec`)
  - `:distance`: the (absolute) distance to the closest point
  - `:signed_distance`: the signed distance to the closest point
- `mesh`: the `SimpleMesh`
- `points`: the point coordinates for which the `fields` are computed; one of the following:
  - `Tuple` with a list of coordinates for each dimension (usually a kind of
      `AbstractRange`)
  - `CartesianGrid`, using the coordinates of the centroids of the grid cells
  - a collection of `Point`s

# Returns

- `NamedTuple` with an `Array` of the same dimension as `points` for each of the `fields`
"""
function rasterize(fields::NTuple{N,Symbol}, mesh::SimpleMesh{Dim,T}, points::Raster{Dim,T},
        method = default_method()) where {N,Dim,T}
    dims = dimensions(points)
    data = map(fields) do field
        if field in (:distance, :signed_distance)
            Array{T,length(dims)}(undef, dims)
        elseif field in (:closest_point, )
            zeros(Point{Dim,T}, dims)
        elseif field in (:direction, )
            zeros(Vec{Dim,T}, dims)
        else
            error("Unsupported field: `$field`")
        end
    end
    if Dim == 3 && !(topology(mesh) isa HalfEdgeTopology)
        mesh = SimpleMesh(vertices(mesh), convert(HalfEdgeTopology, topology(mesh)))
    end
    rasterize!(NamedTuple{fields}(data), mesh, points, method)
end


# -----------------
# Output Management
# -----------------

function current_distance(data, ind)
    if haskey(data, :distance)
        data.distance[ind]
    elseif haskey(data, :signed_distance)
        abs(data.signed_distance[ind])
    else
        error("Data does not contain any field for distance")
    end
end

function reset_distance!(data)
    for key in (:distance, :signed_distance)
        haskey(data, key) || continue
        values = getproperty(data, key)
        fill!(values, typemax(eltype(values)))
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
