module MeshRasterization

include("meshes-piracy.jl")
include("prism.jl")
include("hyperplane.jl")
include("rasters.jl")
include("scan-conversion.jl")
include("rasterization.jl")
include("brute-force.jl")
include("csc.jl")
include("mesh-input.jl")

using .Rasterization, .CSCMethod, .BruteForceMethod, .MeshInput

export closest_point, direction, distance, signed_distance, rasterize, rasterize!,
        BruteForceRasterization, CharacteristicScanConversion, loadmesh

end # module
