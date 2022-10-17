module BruteForceMethod

using LinearAlgebra: norm, dot

# Meshes.jl: basic functionality
using Meshes: Point, Vec, Segment, Triangle, vertices, normal

# Meshes.jl: mesh & topology
using Meshes: SimpleMesh, HalfEdgeTopology, Boundary, Coboundary,
    topology, nfaces, element, facet, indices, connect

import ..Hyperplanes: Hyperplane
import ..Rasterization: RasterizationMethod, rasterize!, pointindices,
    getpoint, current_distance, reset_distance!, vertex_orientation

struct BruteForceRasterization <: RasterizationMethod end

function closest_point(s::Segment{Dim,T}, pt::Point{Dim,T}) where {Dim,T}
    v1, v2 = vertices(s)
    dir = v2 - v1 # direction vector along segment
    l2 = sum(abs2, dir) # length^2 of segment
    loc = clamp(dot(pt - v1, dir) / l2, zero(T), one(T)) # location ∈ (0,1)
    cpt = v1 + loc * dir
    cpt => norm(pt - cpt)
end

function closest_point(hp::Hyperplane{Dim,T}, pt::Point{Dim,T}) where {Dim,T}
    dist = dot(normal(hp), pt - hp.p)
    cpt = pt - dist * normal(hp)
    cpt => abs(dist)
end

hyperplane(t::Triangle) = Hyperplane(first(vertices(t)), normal(t))

function closest_point(t::Triangle{Dim,T}, pt::Point{Dim,T}) where {Dim,T}

    #println("checking point ", pt)

    v1, v2, v3 = vertices(t)

    # first check edge segments
    cpt, dist = closest_point(Segment(Vec(v1, v2)), pt)
    for (s1, s2) in ((v2, v3), (v3, v1))
        cpt_, dist_ = closest_point(Segment(Vec(s1, s2)), pt)
        if dist_ < dist
            #println("  - segment update: ", cpt_)
            dist, cpt = dist_, cpt_
        end
    end

    # then also check distance to plane
    cpt_, dist_ = closest_point(Hyperplane(v1, normal(t)), pt)
    if cpt_ in t
        #println("  - plane update: ", cpt_)
        dist, cpt = dist_, cpt_
    end

    #println("found ", cpt, " at distance ", dist)

    cpt => dist
end

function rasterize!(data::NamedTuple, mesh::SimpleMesh{Dim}, points, ::BruteForceRasterization) where {Dim}
    reset_distance!(data)
    for ind in pointindices(points)
        point = getpoint(points, ind)
        #println("$point:")
        for face in 1:nfaces(mesh, Dim - 1)
            meshpoint, distance = closest_point(element(mesh, face), point)

            # only update if the new distance is shorter
            if distance >= current_distance(data, ind)
                #println(" - skipped $meshpoint at d=$distance")
                continue
            end

            if haskey(data, :closest_point)
                data.closest_point[ind] = meshpoint
            end

            if haskey(data, :direction)
                if meshpoint ≈ point
                    data.direction[ind] = point - point # zero vector
                else
                    dir = meshpoint - point
                    data.direction[ind] = dir / norm(dir)
                end
            end

            if haskey(data, :distance)
                data.distance[ind] = distance
            end

            if haskey(data, :signed_distance)
                T = eltype(data.signed_distance)

                # use `>=` instead of `sign()` to avoid multiplying by zero
                #sign = dot(point - meshpoint, normal(element)) >= zero(T) ? 1 : -1
                type, n = face_orientation(mesh, face, meshpoint)
                sign = dot(point - meshpoint, n) >= zero(T) ? 1 : -1
                #println(" - set o($type)=$sign for $meshpoint at d=$distance")
                data.signed_distance[ind] = distance * sign
            end
        end
    end
    data
end

function within(point::Point{N}, segment::Segment{N}) where {N}
    # or: covers(segment, point), covers(point, point) ?
    closest_point(segment, point)[1] ≈ point
end

function face_orientation(mesh::SimpleMesh{3}, face::Int, point)

    # check if point is located at one of the vertices
    for vertex in Boundary{2,0}(topology(mesh))(face)
        if point ≈ vertices(mesh)[vertex]
            return :vertex, vertex_orientation(mesh, vertex)
        end
    end

    # check if point is located within one of the edges
    for edge in Boundary{2,1}(topology(mesh))(face)
        if within(point, facet(mesh, edge))
            return :edge, edge_orientation(mesh, edge)
        end
    end

    # otherwise return orientation of element itself
    :face, normal(element(mesh, face))
end

function face_orientation(mesh::SimpleMesh{2}, face, point)

    # check if point is located at one of the vertices
    for vertex in element2vertices(mesh, face)
        if point ≈ vertices(mesh)[vertex]
            return :vertex, vertex_orientation(mesh, vertex)
        end
    end

    # otherwise return orientation of element itself
    :face, normal(element(mesh, face))
end

function edge_orientation(mesh::SimpleMesh{3,T,V,<:HalfEdgeTopology}, ind) where {T,V}
    sum(Coboundary{1,2}(topology(mesh))(ind)) do el
        normal(element(mesh, el))
    end
end

function element2vertices(mesh::SimpleMesh{2}, ind)
    indices(element(topology(mesh), ind))
end

function vertex2elements(mesh::SimpleMesh{2}, ind)
    # elements are segments here
    n = length(vertices(mesh))
    e1 = ind == 1 ? (n,1) : (ind-1,ind)
    e2 = ind == n ? (n,1) : (ind,ind+1)
    connect.((e1, e2))
end

end # module
