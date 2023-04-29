module Rasters

using Meshes: CartesianGrid, Point, Vec, spacing as gridspacing, centroid

const RasterAxes{Dim,T} = NTuple{Dim,AbstractRange{T}}
const Raster{Dim,T} = Union{CartesianGrid{Dim,T},RasterAxes{Dim,T}}

# coordinates of point that would have the index zero
origin(raster::CartesianGrid) = minimum(raster) - spacing(raster) / 2
origin(raster::RasterAxes) = Point(first.(raster)) - spacing(raster)

# distance between adjacent points
spacing(raster::CartesianGrid) = Vec(gridspacing(raster))
spacing(raster::RasterAxes) = Vec(step.(raster))

# number of points along each dimension
dimensions(raster::CartesianGrid) = size(raster)
dimensions(raster::RasterAxes) = length.(raster)

# get point at a specific index
getpoint(raster::CartesianGrid, ind) = centroid(raster, ind)
getpoint(raster::RasterAxes, ind) = Point(getindex.(raster, Tuple(ind)))

#= other raster functions, currently unused
getpoint(points, ind) = points[ind]
pointindices(points) = eachindex(points)
pointindices(points::Tuple) = CartesianIndices(eachindex.(points))
pointindices(grid::CartesianGrid) = CartesianIndices(size(grid))
dimensionality(raster::CartesianGrid{Dim}) where Dim = Dim
dimensionality(raster::NTuple{N,AbstractRange{T}}) where {N,T} = N
# =#

end # module Rasters
