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


for raster in ("CartesianGrid", "ranges")
    @testset "$(string(method)) for $(raster)" for method in methods
        @testset "$(case.label)" for case in cases2d
            ng = Ngon(case.points...)
            gd = if raster == "CartesianGrid"
                CartesianGrid(case.grid[1], convert.(Float64, case.grid[2]),
                              convert.(Float64, case.grid[3]), case.grid[4:end]...)
            elseif raster == "ranges"
                ox, oy = length(case.grid) > 3 ? case.grid[4] : (1, 1)
                lx, ly = case.grid[1]
                dx, dy = case.grid[3]
                x0, y0 = case.grid[2] .+ (1 .- (ox, oy)) .* (dx, dy)
                rx = range(x0+dx/2, x0+lx-dx/2, round(Int, lx/dx))
                ry = range(y0+dy/2, y0+ly-dy/2, round(Int, ly/dy))
                (rx, ry)
            end
            broken = method in get(case, :broken, ())
            @test [x.I for x in scan(ng, gd, method())] ==
            collect(case.result) broken=broken
        end
    end
end
