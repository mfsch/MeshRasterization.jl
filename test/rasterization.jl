const fields = (:closest_point, :direction, :distance, :signed_distance)

@testset "Find closest point" begin
    points = [(1,1),(2,1),(1,2)]
    connec = [(1,2),(2,3),(3,1)]
    mesh = SimpleMesh(points, connect.(connec))
    pt = Point(1.25, 1.75)
    r = [MeshRasterization.BruteForceMethod.closest_point(s, pt) for s in mesh]
    @test r[1][1] ≈ Point(1.25, 1)
    @test r[2][1] ≈ Point(1.25, 1.75)
    @test r[3][1] ≈ Point(1, 1.75)
    @test r[1][2] ≈ 0.75
    @test r[2][2] + 1 ≈ 1
    @test r[3][2] ≈ 0.25
end

@testset "Simple 2D rasterization" begin

    points = [(1,1),(2,1),(1,2)]
    connec = [(1,2),(2,3),(3,1)]
    mesh = SimpleMesh(points, connect.(connec))
    grid = CartesianGrid((0.0,0.0), (2.0,2.0), dims=(4,4))
    r = rasterize(fields, mesh, grid)

    @test all(haskey(r, f) for f in fields)

    # check closest points
    @test all(r.closest_point[i,j] ≈ Point(1,1) for i=1:2, j=1:2)
    @test all(r.closest_point[i,3] ≈ Point(1,1.25) for i=1:2)
    @test all(r.closest_point[i,4] ≈ Point(1,1.75) for i=1:2)
    @test all(r.closest_point[3,j] ≈ Point(1.25,1) for j=1:2)
    @test all(r.closest_point[4,j] ≈ Point(1.75,1) for j=1:2)
    @test r.closest_point[3,3] ≈ Point(1,1.25) || r.closest_point[3,3] ≈ Point(1.25,1)
    @test r.closest_point[3,4] ≈ Point(1.25,1.75)
    @test r.closest_point[4,3] ≈ Point(1.75,1.25)
    @test r.closest_point[4,4] ≈ Point(1.5,1.5)

    # check distances
    d1 = sqrt(0.25^2 + 0.25^2)
    d2 = sqrt(0.25^2 + 0.75^2)
    d3 = sqrt(0.75^2 + 0.75^2)
    @test r.distance ≈ [d3 d2 0.75 0.75; d2 d1 0.25 0.25;
                        0.75 0.25 0.25 0.0; 0.75 0.25 0.0 d1]

    # check signed distances
    for i=1:4, j=1:4
        sign = (i == j == 3) ? -1 : 1
        # add 1 to avoid comparing zeros
        @test r.distance[i,j] * sign + 1 ≈ r.signed_distance[i,j] + 1
    end

    # check direction vectors
    d(x,y) = (n = sqrt(x^2 + y^2); Vec(x/n, y/n)) # helper function
    @test r.direction[1:2,1:2] ≈ hcat([d(1,1), d(1,3)], [d(3,1), d(1,1)])
    @test r.direction[3:4,1:2] ≈ ones(2,2) .* (d(0,1),)
    @test r.direction[1:2,3:4] ≈ ones(2,2) .* (d(1,0),)
    @test r.direction[3,3] ≈ d(-1,0) || r.direction[3,3] ≈ d(0,-1)
    @test r.direction[3,4] ≈ Vec(0,0)
    @test r.direction[4,3] ≈ Vec(0,0)
    @test r.direction[4,4] ≈ d(-1,-1)
end


@testset "Degenerate 2D rasterization" begin
    points = [(0,0), (2,0), (2,1), (3,0), (5,0), (5,5), (3,5), (2,4), (2,5), (0,5)]
    connec = [i == length(points) ? (i,1) : (i,i+1) for i=1:length(points)]
    mesh = SimpleMesh(points, connect.(connec))
    grid = (1:4.0, 1:4.0)
    r = rasterize(fields, mesh, grid)
    @test r.signed_distance ≈ - r.distance
end

@testset "Simple 3D rasterization" begin

    points = [(1,1,1),(2,1,1),(1,2,1),(1,1,2)]
    connec = [(3,2,1),(4,3,1),(1,2,4),(2,3,4)]
    mesh = SimpleMesh(points, connect.(connec))
    grid = CartesianGrid((0.,0.,0.), (2.,2.,2.), dims=(4,4,4))
    r = rasterize(fields, mesh, grid)

    for k=1:2
        for i=1:2, j=1:2
            @test r.closest_point[i,j,k] ≈ Point(points[1])
        end
        for i=1:2
            @test r.closest_point[i,3,k] ≈ Point(1, 1.25, 1)
        end
        for j=1:2
            @test r.closest_point[3,j,k] ≈ Point(1.25, 1, 1)
        end
        @test r.closest_point[3,3,k] ≈ Point(1.25, 1.25, 1)
        @test r.closest_point[3,4,k] ≈ Point(1.25, 1.75, 1)
        @test r.closest_point[4,3,k] ≈ Point(1.75, 1.25, 1)
        @test r.closest_point[4,4,k] ≈ Point(1.5, 1.5, 1)
    end
    for i=1:2, j=1:2
        @test r.closest_point[i,j,3] ≈ Point(1, 1, 1.25)
        @test r.closest_point[i,j,4] ≈ Point(1, 1, 1.75)
    end
    for i=1:2
        @test r.closest_point[i,3,3] ≈ Point(1, 1.25, 1.25)
        @test r.closest_point[i,4,3] ≈ Point(1, 1.75, 1.25)
    end
    for j=1:2
        @test r.closest_point[3,j,3] ≈ Point(1.25, 1, 1.25)
        @test r.closest_point[4,j,3] ≈ Point(1.75, 1, 1.25)
    end

    @test r.distance[1,1,1] ≈ sqrt(3) * 0.75

    for i=1:4, j=1:4, k=1:4
        sign = (i == j == k == 3) ? -1 : 1
        @testset "i=$i, j=$j, k=$k" begin
            @test sign * r.distance[i,j,k] ≈ r.signed_distance[i,j,k]
        end
    end
end

@testset "Larger 3D rasterization" begin

    points = [(0,0,0),(3,1,1),(1,3,1),(1,1,3)]
    connec = [(3,2,1),(4,3,1),(1,2,4),(2,3,4)]
    mesh = SimpleMesh(points, connect.(connec))

    N = 10
    grid = CartesianGrid((0.,0.,0.), (3.,3.,3.), dims=(N, N, N))
    r = rasterize(fields, mesh, grid)

    dx = 3 / N
    for i=2:N-1, j=2:N-1, k=2:N-1
        if r.distance[i,j,k] > sqrt(3) * dx
            outside = r.signed_distance[i,j,k] >= 0
            @test all(outside == (r.signed_distance[i+di,j+dj,k+dk] >= 0)
                      for di=-1:1, dj=-1:1, dk=-1:1)
        end
    end

end

@testset "CSC: 2D rasterization" begin #=
    points = [(1,1),(2,1),(1,2)]
    connec = [(1,2),(2,3),(3,1)]
    mesh = SimpleMesh(points, connect.(connec))
    grid = CartesianGrid((0.0,0.0), (2.0,2.0), dims=(4,4))

    rb = rasterize(fields, mesh, grid)
    rc = rasterize(fields, mesh, grid, CharacteristicScanConversion(2.0))

    @test rb.distance ≈ rc.distance
    @test rb.signed_distance ≈ rc.signed_distance
    @test rb.direction ≈ rc.direction
    @test rb.closest_point == rc.closest_point # points are compared approximately
=# end

@testset "CSC: 3D rasterization" begin
    points = [(1,1,1),(2,1,1),(1,2,1),(1,1,2)]
    connec = [(3,2,1),(4,3,1),(1,2,4),(2,3,4)]
    mesh = SimpleMesh(points, connect.(connec))
    grid = CartesianGrid((0.,0.,0.), (2.,2.,2.), dims=(4,4,4))

    rb = rasterize(fields, mesh, grid)
    rc = rasterize(fields, mesh, grid, CharacteristicScanConversion(2.0))

    @test rb.distance ≈ rc.distance
    @test rb.signed_distance ≈ rc.signed_distance
    @test rb.direction ≈ rc.direction
    @test rb.closest_point == rc.closest_point # points are compared approximately
end

@testset "Coordinates on boundary" begin
    pts = [(0, 1, 0), (0.5, 0.5, 0), (0, 0, 0), (1, 0, 0), (1, 1, 0)]
    els = [(1, 2, 3), (2, 4, 3), (1, 5, 4)]
    mesh = SimpleMesh([Float64.(pt) for pt in pts], connect.(els))
    x = range(0, 1, 32+1)[1:end-1]
    raster = (x, x, 0:1.0:0)
    rb = signed_distance(mesh, raster)
    rc = signed_distance(mesh, raster, CharacteristicScanConversion(1))
    @test rb ≈ rc
end

@testset "Performance checks" begin
    points = [(1,1,1),(2,1,1),(1,2,1),(1,1,2)]
    connec = [(3,2,1),(4,3,1),(1,2,4),(2,3,4)]
    mesh = SimpleMesh(points, connect.(connec), relations=true)
    grid = CartesianGrid((0.,0.,0.), (2.,2.,2.), dims=(100, 100, 100))

    dims = size(grid)
    data = (distance = zeros(dims),
            signed_distance = zeros(dims),
            closest_point = zeros(Point{3,Float64}, dims),
            direction = zeros(Vec{3,Float64}, dims))

    fill!(data.distance, Inf)
    tb = rasterize!(data, mesh, grid, BruteForceRasterization())
    fill!(data.distance, Inf)
    tc = rasterize!(data, mesh, grid, CharacteristicScanConversion(2.0))
    fill!(data.distance, Inf)
    tb = @elapsed rasterize!(data, mesh, grid, BruteForceRasterization())
    fill!(data.distance, Inf)
    tc = @elapsed rasterize!(data, mesh, grid, CharacteristicScanConversion(2.0))
    println("time: $tb for brute force, $tc for CSC")
    @test tb > tc
end
