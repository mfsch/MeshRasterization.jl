module CSCMethod

using LinearAlgebra: norm, dot, cross
using Meshes: Vec, Ball, vertices, normal, isconvex, SimpleMesh, Boundary,
              Coboundary, topology, nfaces, element
using ..Hyperplanes: Hyperplane
using ..Prisms: Prism
using ..Rasters: Raster, getpoint
using ..ScanConversion: EdgeFunctionScan, scan
using ..Rasterization: RasterizationMethod, current_distance, reset_distance!,
                       vertex_orientation

import ..Rasterization: rasterize!

export CharacteristicScanConversion

struct CharacteristicScanConversion{T} <: RasterizationMethod
    dmax::T
end

function rasterize!(data::NamedTuple, mesh::SimpleMesh{Dim,T}, points::Raster{Dim,T},
        method::CharacteristicScanConversion; verbose = false) where {Dim,T}
    reset_distance!(data)

    dmax = abs(method.dmax)
    topo = topology(mesh)

    # plane extrusions
    for face in 1:nfaces(topo, 2)
        ngon = element(mesh, face)
        @assert isconvex(ngon)
        n = normal(ngon)
        v = vertices(ngon)
        prism = Prism(v .- (dmax * n,), 2 * dmax * n)

        for ind in scan(prism, points, EdgeFunctionScan())
            point = getpoint(points, ind)
            signed_distance = dot(point - first(v), n)
            distance = abs(signed_distance)

            # only update if the new distance is shorter
            if distance > dmax || distance >= current_distance(data, ind)
                continue
            end

            if haskey(data, :closest_point)
                data.closest_point[ind] = point - n * signed_distance
            end

            if haskey(data, :direction)
                if distance + 1 ≈ 1
                    data.direction[ind] = zero(Vec{Dim,T})
                else
                    data.direction[ind] = - sign(signed_distance) * n
                end
            end

            if haskey(data, :distance)
                data.distance[ind] = distance
            end

            if haskey(data, :signed_distance)
                data.signed_distance[ind] = signed_distance
            end
        end
    end

    # edge extrusions
    edge2faces = Coboundary{1,2}(topo)
    edge2vertices = Boundary{1,0}(topo)
    for edge in 1:nfaces(topo, 1)
        planar_edges = 0 # only for monitoring

        faces = edge2faces(edge)
        if length(faces) == 1
            # TODO: handle border edges
            verbose && println("skipped border edge")
            continue
        end

        n1, n2 = normal.(element.((mesh,), faces)) # should always have 2 neighbors

        if n1 ≈ n2
            planar_edges += 1
            continue # planar edges are covered by plane extusions
        end

        v1, v2 = vertices(mesh)[edge2vertices(edge)]
        sign = dot(cross(v2 - v1, n1), n2) > zero(T) ? 1 : -1 # positive for convex edge

        # compute base length of prism such that the middle is `dmax` away from the edge
        # TODO: use multiple prisms for edges with acute angles
        ratio = norm(n1 + n2) / (1 + dot(n1, n2)) # 1 / cos(α)
        ratio < 10 * one(T) || error("Degenerate mesh: acute edge angle")
        lbase = dmax * ratio * sign

        # define prism for edge extrusion
        base = Vec(v1, v1 + lbase * n1, v1 + lbase * n2)
        side = v2 - v1
        prism = Prism(base, side)

        edgedir = side / norm(side)
        for ind in scan(prism, points, EdgeFunctionScan())
            point = getpoint(points, ind)
            distance = norm(cross(point - v1, edgedir))
            signed_distance = sign * distance

            # only update if the new distance is shorter
            if distance > dmax || distance >= current_distance(data, ind)
                continue
            end

            meshpoint = if haskey(data, :closest_point) || haskey(data, :direction)
                v1 + dot(point - v1, edgedir) * edgedir
            else
                nothing # do not compute if not needed
            end

            if haskey(data, :closest_point)
                data.closest_point[ind] = meshpoint
            end

            if haskey(data, :direction)
                if distance + 1 ≈ 1
                    data.direction[ind] = zero(Vec{Dim,T})
                else
                    dir = meshpoint - point
                    data.direction[ind] = dir / norm(dir)
                end
            end

            if haskey(data, :distance)
                data.distance[ind] = distance
            end

            if haskey(data, :signed_distance)
                data.signed_distance[ind] = signed_distance
            end
        end

    end
    verbose && println("skipped $border_edges border edges, $flat_edges flat edges")

    # vertex extrusions
    for vertex in 1:nfaces(topo, 0)

        meshpoint = vertices(mesh)[vertex]

        # define ball for vertex extrusion
        ball = Ball(meshpoint, dmax)

        for ind in scan(ball, points)
            point = getpoint(points, ind)
            dir = meshpoint - point # not normalized!
            distance = norm(dir)

            # only update if the new distance is shorter
            if distance > dmax || distance >= current_distance(data, ind)
                continue
            end

            if haskey(data, :closest_point)
                data.closest_point[ind] = meshpoint
            end

            if haskey(data, :direction)
                if distance + 1 ≈ 1
                    data.direction[ind] = zero(Vec{Dim,T})
                else
                    data.direction[ind] = dir / distance
                end
            end

            if haskey(data, :distance)
                data.distance[ind] = distance
            end

            if haskey(data, :signed_distance)
                # sign is flipped because dir is oriented from point to vertex
                sgn = dot(dir, vertex_orientation(mesh, vertex)) < 0 ? 1 : -1
                data.signed_distance[ind] = sgn * distance
            end
        end
    end

    data
end

end # module
