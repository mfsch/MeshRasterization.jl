using Test, Meshes, SDF

@testset "SDF.jl" begin
    @testset "Scan Conversion" begin include("scan-conversion.jl") end
    @testset "Rasterization" begin include("rasterization.jl") end
end
