module MeshRasterization

include("meshes-piracy.jl")
include("prism.jl")
include("hyperplane.jl")
include("scan-conversion.jl")
include("rasterization.jl")
include("brute-force.jl")
include("csc.jl")

using Meshes: Point, Vec, SimpleMesh, HalfEdgeTopology, vertices, topology
import .Rasterization: rasterize, rasterize!, pointdims
import .BruteForceMethod: BruteForceRasterization
import .CSCMethod: CharacteristicScanConversion
export closest_point, direction, distance, signed_distance, rasterize, rasterize!,
        BruteForceRasterization, CharacteristicScanConversion


# --------------------
# High-level Interface
# --------------------

function closest_point(mesh, points)
    rasterize((:closest_point,), mesh, points).closest_point
end

function direction(mesh, points)
    rasterize((:direction,), mesh, points).direction
end

function distance(mesh, points)
    rasterize((:distance,), mesh, points).distance
end

function signed_distance(mesh, points)
    rasterize((:signed_distance,), mesh, points).signed_distance
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
function rasterize(fields::NTuple{N,Symbol}, mesh::SimpleMesh{Dim,T}, points,
        method = BruteForceRasterization()) where {N,Dim,T}
    dims = pointdims(points)
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

end # module
