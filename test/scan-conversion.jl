using MeshRasterization.ScanConversion

const methods = (SimpleScan, EdgeFunctionScan)
const cases2d = [ # grid = (dims, origin, spacing[, offset])
                    (label = "unit square",
                     grid = ((2,2), (-1,-1), (1,1)),
                     points = ((0,0), (1,0), (1,1), (0,1)),
                     result = ((2,2),)),
                    (label = "unit square with offset",
                     grid = ((2,2), (-1,-1), (1,1), (0,0)),
                     points = ((0,0), (1,0), (1,1), (0,1)),
                     result = ((1,1),)),
                    (label = "simple triangle",
                     grid = ((3,4), (0,0), (1/3,1/4)),
                     points = ((1,0), (1,1), (0,1)),
                     result = ((3,2), (2,3), (3,3), (1,4), (2,4), (3,4))),
                    (label = "triangle with border points",
                     grid = ((3,3), (0,0), (1/3,1/3)),
                     points = ((1,0), (1,1), (0,1)),
                     result = ((3,1), (2,2), (3,2), (1,3), (2,3), (3,3)),
                     broken = ()), # TraverseBoundingBox may or may not pass
                    (label = "small triangle",
                     grid = ((3,4), (0,0), (1,1)),
                     points = ((1.25,2.25), (1.75,2.25), (1.5,2.75)),
                     result = ((2,3),)),
                   ]



@testset "$(string(method))" for method in methods
    @testset "$(case.label)" for case in cases2d
        ng = Ngon(case.points...)
        gd = CartesianGrid(case.grid[1],
                           convert.(Float64, case.grid[2]),
                           convert.(Float64, case.grid[3]), case.grid[4:end]...)
        broken = method in get(case, :broken, ())
        @test [x.I for x in scan(ng, gd, method())] ==
        collect(case.result) broken=broken
    end
end
