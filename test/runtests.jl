using Test, Meshes, MeshRasterization

@testset "MeshRasterization.jl" begin
    @testset "Scan Conversion" begin include("scan-conversion.jl") end
    @testset "Rasterization" begin include("rasterization.jl") end
end
