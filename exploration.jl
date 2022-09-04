### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 9929f2be-277a-11ed-375a-f18a1b52c0ca
begin
	LOAD_PATH = ["@", "@stdlib"]
	import Pkg
	Pkg.activate(".")
	import LinearAlgebra, FileIO, MeshViz, Meshes, JSServe, Profile, PProf
	import WGLMakie as Mke
	using Test, BenchmarkTools
	JSServe.configure_server!(listen_url = "100.85.63.104", listen_port = 8146, external_url = "http://100.85.63.104:8146")
	isassigned(JSServe.GLOBAL_SERVER) && close(JSServe.GLOBAL_SERVER[])
end;

# ╔═╡ 87aca29f-29df-44ba-9630-cc97e927a339
present(x) = (show(stdout, MIME"text/plain"(), x); println())

# ╔═╡ 3671d217-c327-431f-9bfd-0f3e7fd4a01a
md"# Generate SDF from STL files"

# ╔═╡ fc1826d8-f6ad-443f-bf65-2d5def8974a7
md"## Load STL from file"

# ╔═╡ c4697b84-9cbd-4178-9aaf-72527bddde67
function loadmesh(path) # from MeshBridge.jl
	mesh = FileIO.load(path)

	vertexToIdx = Dict()
    for i in 1:length(mesh.position)
        vertex = mesh.position[i]
        if haskey(vertexToIdx, vertex)
            continue
        end
        vertexToIdx[vertex] = i
    end
    faces = []
    for triangle in mesh
        i1 = vertexToIdx[triangle[1]]
        i2 = vertexToIdx[triangle[2]]
        i3 = vertexToIdx[triangle[3]]
        push!(faces, (i1, i2, i3))
    end
    topology = Meshes.FullTopology(Meshes.connect.(faces))
    result = Meshes.SimpleMesh(
        [Meshes.Point([x[1], x[2], x[3]]) for x in mesh.position], topology
    )
    return result
end

# ╔═╡ f33325a0-ca9e-4366-bc37-91dc414cad25
let mesh = loadmesh("single-cube.stl")
	MeshViz.viz(mesh, showfacets=true, color=1:Meshes.nelements(mesh))
end

# ╔═╡ 34b3f98c-1053-4aa7-ae0f-4aff49998c82
let dims = (3, 3, 3), origin = (0.25, 0.25, 0.25), spacing = (0.5, 0.5, 0.5)
	gd = Meshes.CartesianGrid(dims, origin, spacing)
	gd = Meshes.vertices(gd)
	MeshViz.viz(gd, size=10)
end

# ╔═╡ 82d773d2-a328-41ca-b5c6-1836a393d6da
md"## Rasterize 2D Polygons"

# ╔═╡ a45d21c3-ea3b-4d11-8a4a-be519f1d7719
begin
	abstract type Rasterization end
	struct TraverseGrid <: Rasterization end
	struct TraverseBoundingBox <: Rasterization end
	struct EdgeFunctions{T} <: Rasterization
		tol::T
		EdgeFunctions(tol = 1e-9) = new{typeof(tol)}(tol)
	end

	# test if point in polygon for all points in grid
	function rasterize(ngon::Meshes.Ngon{N,2,T},
					   grid::Meshes.CartesianGrid{2,T},
					   ::TraverseGrid;
					   verbose = false) where {N,T}
		verbose && println("Rasterize by traversing $(join(size(grid), '×')) grid:")
		
		lininds = LinearIndices(size(grid))
		Iterators.filter(CartesianIndices(size(grid))) do i
			pt = Meshes.centroid(grid[lininds[i]])
			if pt ∈ ngon
				verbose && println(" ☒ $(i.I): $pt")
				true
			else
				verbose && println(" ☐ $(i.I): $pt")
				false
			end
		end
	end

	function rasterize(box::Meshes.Box{N,T}, grid::Meshes.CartesianGrid{N,T}) where {N,T}
		dx = Meshes.spacing(grid)
		imin = (minimum(box) - minimum(grid)) ./ dx .+ 1|>
			i -> floor.(Int, i) |>
			i -> max.(i, (1, 1))
		imax = (maximum(box) - minimum(grid)) ./ dx .+ 1|>
			i -> floor.(Int, i) |>
			i -> min.(i, size(grid))
		inds = ntuple(i -> imin[i]:imax[i], N)
		CartesianIndices(inds)
	end

	# test if point in polygon for all points in bounding box
	function rasterize(ngon::Meshes.Ngon{N,2,T},
					   grid::Meshes.CartesianGrid{2,T},
		     		   ::TraverseBoundingBox;
					   verbose = false) where {N,T}

		Δgd = Meshes.spacing(grid)
		xref = minimum(grid) - 0.5 * Δgd

		inds = rasterize(Meshes.boundingbox(ngon), grid)
		verbose && println("Rasterize by traversing $(join(size(inds), '×')) bounding box:")
		
		Iterators.filter(inds) do i
			pt = xref + Δgd .* i.I
			if pt ∈ ngon
				verbose && println(" ☒ $(i.I): $pt")
				true
			else
				verbose && println(" ☐ $(i.I): $pt")
				false
			end
		end
	end

	# test points using edge functions (Pineda, 1988)
	# → only works for convex polygons!
	function rasterize(ngon::Meshes.Ngon{N,2,T},
					   grid::Meshes.CartesianGrid{2,T},
		     		   method::EdgeFunctions;
					   verbose = false) where {N,T}
		
		Δgd = Meshes.spacing(grid)
		X = Meshes.vertices(ngon)
		dX = ntuple(i -> X[i] - X[i == 1 ? N : i-1], N)

		# Ei for reference position: centroid of index (0,0)
		xref = minimum(grid) - 0.5 * Δgd
		Eref = ntuple(i -> Meshes.cross(xref-X[i], dX[i]), N)

		# change of Ei for moving one point in x-direction
		dEx = ntuple(i -> dX[i][2] * Δgd[1], N)

		# change of Ei for moving one point in y-direction
		dEy = ntuple(i -> - dX[i][1] * Δgd[2], N)

		# check edge functions for all indices in bounding box
		inds = rasterize(Meshes.boundingbox(ngon), grid)
		check_edgefunctions(Eref, dEx, dEy, inds; tol = convert(T, method.tol))
	end

	function check_edgefunctions(Eref::NTuple{N,T}, dEx::NTuple{N,T}, dEy::NTuple{N,T}, inds::CartesianIndices{2}; tol::T = eps(T)) where {N,T}
		Iterators.filter(inds) do ind
			all(Eref[i] + ind[1] * dEx[i] + ind[2] * dEy[i] < tol for i in 1:N)
		end
	end

	rasterize
end

# ╔═╡ 07d89b37-ffc9-49c1-9f7c-59e26f256fcd
let
	#pg = Meshes.Ngon([(1,8+1e-9), (5+1*rand(),8), (4,2), (9.5-1e-9,5+1e-9), (1,1.0)])

	# set up polygon
	pts = Meshes.PointSet([(1,8+1e-9), (5+1*rand(),8), (4,2), (9.5-1e-9,5+1e-9), (1,1.0)])
	pg = Meshes.hull(pts, Meshes.GrahamScan()) |> Meshes.vertices |> collect |> Meshes.Ngon

	N = (32, 32)

	
	gd = Meshes.CartesianGrid(N, Meshes.Point(0, 0), Meshes.Vec(10 ./ N))
	data = Meshes.meshdata(gd, Dict(2 => (distance=Union{Missing, Float64}[missing for el in gd], )))
	

	for ind in rasterize(pg, gd, EdgeFunctions())
		i = LinearIndices(size(gd))[ind]
		data["distance"][i] = 0.5
	end
	
	fig = Mke.Figure(resolution = (400, 400))
	ax = Mke.Axis(fig[1,1])

	# add distance values to figure
	col = map(data["distance"]) do d
		ismissing(d) ? "rgb(240, 240, 240)" : MeshViz.colorschemes[:viridis][d]
	end
	MeshViz.viz!(ax, gd, color=col, showfacets=true)
	MeshViz.viz!(ax, Meshes.centroid.(gd), color=:white, size=5)

	# add polygon to figure
	MeshViz.viz!(ax, pg, alpha=0.5)
	MeshViz.viz!(ax, pts, color=:black)

	fig
end

# ╔═╡ 17854bbc-543c-42bd-a126-943edda8c75a
md"### Tests"

# ╔═╡ 3772c9e9-2db1-4978-a011-de19799f4e98
const test_cases = [ # grid = (dims, origin, spacing[, offset])
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
	 broken = (TraverseGrid, )), # TraverseBoundingBox may or may not pass
	(label = "small triangle",
	 grid = ((3,4), (0,0), (1,1)),
	 points = ((1.25,2.25), (1.75,2.25), (1.5,2.75)),
	 result = ((2,3),)),
]

# ╔═╡ cb2c81fe-1987-42c9-aec8-badbbd0b1a0a
function runtests(cases, methods)
	@testset "$(split(string(method), '.')[end])" for method in methods
		@testset "$(case.label)" for case in cases
			ng = Meshes.Ngon(case.points...)
			gd = Meshes.CartesianGrid(case.grid[1],
				convert.(Float64, case.grid[2]),
				convert.(Float64, case.grid[3]), case.grid[4:end]...)
			broken = method in get(case, :broken, ())
			@test [x.I for x in rasterize(ng, gd, method(), verbose=false)] ==
				collect(case.result) broken=broken
		end
	end
end

# ╔═╡ 2d4b6ea3-decb-4f6a-939d-1eafe956828e
@testset "Rasterization" begin
	runtests(test_cases, (TraverseGrid, TraverseBoundingBox, EdgeFunctions))
end;

# ╔═╡ 0d7b073f-90a3-4588-b4d1-5cc0b14e2bc5
md"### Benchmarks"

# ╔═╡ 409498df-51ed-49e8-85a9-efcdee401255
function runbenchmarks(method)
	ngon = Meshes.Ngon([(1.2,1.7), (1.8,1.8), (1.3,2.1), (1.1,1.9)])
	grid = Meshes.CartesianGrid((256,256), (0.,0.), (4π/256, 2π/256))
	@benchmark maximum(rasterize($ngon, $grid, $method))
	#rasterize(ngon, grid, method)
end

# ╔═╡ ef3d6f6e-327f-45b8-92c6-423dd52bac49
# ╠═╡ disabled = true
#=╠═╡
runbenchmarks(TraverseGrid())
  ╠═╡ =#

# ╔═╡ e0415fec-1d32-458d-ae55-87251486cbcb
# ╠═╡ disabled = true
#=╠═╡
runbenchmarks(TraverseBoundingBox())
  ╠═╡ =#

# ╔═╡ 01843255-f9a3-49a7-a989-b15b2e6932de
# ╠═╡ disabled = true
#=╠═╡
runbenchmarks(EdgeFunctions())
  ╠═╡ =#

# ╔═╡ 3c555464-59e2-4c98-bdfe-f64891f6b569
md"### Allocations"

# ╔═╡ bf425675-9ee5-43de-a67a-b9098a751187
# ╠═╡ disabled = true
#=╠═╡
let
	ngon = Meshes.Ngon([(1.2,1.7), (1.8,1.8), (1.3,2.1), (1.1,1.9)])
	grid = Meshes.CartesianGrid((256,256), (0.,0.), (4π/256, 2π/256))
	#method = TraverseBoundingBox()
	method = EdgeFunctions()
	maximum(rasterize(ngon, grid, method), init = CartesianIndex(1,1))
	Profile.Allocs.clear()
	Profile.Allocs.start(sample_rate=1.0)
	maximum(rasterize(ngon, grid, method), init = CartesianIndex(1,1))
	Profile.Allocs.stop()
	PProf.Allocs.pprof(web=false)
end
  ╠═╡ =#

# ╔═╡ 342ed4c3-ffcc-4a75-887d-5c9c6f8318d6
md"## Rasterize 3D Polyhedron"

# ╔═╡ 73dfeea2-3faf-458c-ad41-7d592a581479


# ╔═╡ 1f405903-9e56-408c-aea6-1f1328bc58c9


# ╔═╡ 365fd982-261a-4c11-85d3-0f54b6d4865a


# ╔═╡ 3d03fcd5-a804-452e-b409-8d7317df014f


# ╔═╡ 1b06df66-f0db-40c6-a8ad-79e855c883e0


# ╔═╡ 3bd42cac-c944-44c3-a8d9-e4eb07d059d3


# ╔═╡ a1588120-35c9-4719-8bca-757f56cc288a


# ╔═╡ ff1d3625-60e0-41dd-b889-0d2606378568
md"""# Old Code"""

# ╔═╡ a07eeeed-1d9d-4574-bad3-fabc4f0b0af9
md"""## Find coordinates inside 2D polygon"""

# ╔═╡ f809754a-ea9f-4062-8a94-85e64b42c3a2
begin
	# pts is a selection of points (in any order) defining a convex polygon
	# returns the integer values inside the polygon
	# to find indices, convert the polygon points into a space of continuous indices
	function scan_polygon(pts::Array{Tuple{Float64,Float64},1}, ε = 1e-6)
	    
	    # list of integer coordinates that are contained in the polygon
	    coords = Tuple{Int,Int}[]
	    
	    # find y-extent of polygon
	    sort!(pts, by=pt->(pt[2],pt[1])) # sort by y, then by x
	    
	    # line from imin to imax, in the form v0 + vl·t
	    v0 = [pts[1]...] # vector to origin of the line
	    vl = [pts[end]...] - v0 # vector along the line
	    vl /= (n=LinearAlgebra.norm(vl); n>ε ? n : 1)
	    dx(pt) = pt[1] - v0[1] - LinearAlgebra.dot([pt...]-v0, vl) * vl[1]
	    
	    # split points as below (-1), on (0) or above (1) the line connecting first & last point
	    xdir = map(pt -> (dx_pt=dx(pt); abs(dx_pt)<ε ? 0 : sign(dx_pt)), pts)
	    pts_lower = pts[xdir .<= 0]
	    pts_upper = pts[xdir .>= 0]
	    
	    il, iu = 2, 2 # set indices to second node (end of first line)
	    
	    for y = ceil(Int, pts[1][2]-ε):floor(Int, pts[end][2]+ε)
	        
	        # update il & iu
	        while (pts_lower[il][2] < y && il < length(pts_lower))
	            il += 1
	        end
	        while (pts_upper[iu][2] < y && iu < length(pts_upper))
	            iu += 1
	        end
	        
	        # find extent of x-indices
	        xmin = ((x1,y1) = pts_lower[il-1]; (x2,y2) = pts_lower[il]; 
	            abs(y2-y1)<ε ? min(x1, x2) : x1 + (x2-x1) * (y-y1)/(y2-y1))
	        xmax = ((x1,y1) = pts_upper[iu-1]; (x2,y2) = pts_upper[iu]; 
	            abs(y2-y1)<ε ? max(x1, x2) : x1 + (x2-x1) * (y-y1)/(y2-y1))
	        
	        for x = ceil(Int, xmin-ε):floor(Int, xmax+ε)
	            push!(coords, (x,y))
	        end
	    end
	    
	    coords
	end
	
	scan_polygon(pts::Array{Tuple{T,T},1}) where T <: Real = scan_polygon([convert(Tuple{Float64,Float64}, pt) for pt=pts])

end

# ╔═╡ d4065db9-6b48-4c1e-ae33-ed793fc91c28
@testset "Degenerate Cases 2D" begin
	# test some degenerate cases
	@test all(scan_polygon([(0,0), (0,0), (0,0)]) .== [(0,0)])
	@test all(scan_polygon([(0,0), (0,2)]) .== [(0,0), (0,1), (0,2)])
	@test all(scan_polygon([(0,0), (2,0)]) .== [(0,0), (1,0), (2,0)])
	@test all(scan_polygon([(0,0), (1,1), (1,0), (0,1)]) .== [(0,0), (1,0), (0,1), (1,1)])
	@test (ε=1e-9; all(scan_polygon([(0+ε,0+ε), (1-ε,1-ε), (1-ε,0+ε), (0+ε,1-ε)]) .== [(0,0), (1,0), (0,1), (1,1)]))
end

# ╔═╡ f415eab0-b090-495c-a79c-64a8ad4f3232
function test_scan_polygon(pts)

    # prepare figure
    #figure(figsize=(6,6))
	fig = Mke.Figure(resolution=(600, 600))
	ax = Mke.Axis(fig[1, 1])
    #xlim(-0.25,10.25)
    #ylim(-0.25,10.25)

    # compute output (in ix-iy space) & add to data
    indices = scan_polygon([(1 + 2*pt[1], 1 + 2*pt[2]) for pt=pts])

    # plot data on pixel grid
    data = [0 for i=0:20, j=0:20]
    for i in indices
        if (1 <= i[1] <= size(data,1) && 1 <= i[2] <= size(data,2))
            data[i...] = 1
        else
            println("skipped invalid index: ", i)
        end
    end
    xplt = range(-0.25,10.25,22)
    yplt = range(-0.25,10.25,22)
    #ax[:pcolormesh](xplt, yplt, data')
	Mke.heatmap!(ax, xplt, yplt, data)

    # draw pixel grid midpoints
    xvals = (xplt[1:end-1]+xplt[2:end])/2
    yvals = (yplt[1:end-1]+yplt[2:end])/2
	
    #ax[:plot]([x for x=xvals, y=yvals][:], [y for x=xvals, y=yvals][:], ".w", ms=1)
    Mke.scatter!(ax, [x for x=xvals, y=yvals][:], [y for x=xvals, y=yvals][:], color=:white, markersize=3)
	
    # draw polygon points
    #ax[:plot]([pt[1] for pt=pts], [pt[2] for pt=pts], ".k", ms=10)
    Mke.scatter!(ax, [pt[1] for pt=pts], [pt[2] for pt=pts], color=:black, markersize=10)
    Mke.lines!(ax, [pt[1] for pt=pts], [pt[2] for pt=pts], color=:black)

	fig
end

# ╔═╡ 762e6049-bf06-437c-b990-590a2d025021
test_scan_polygon([(1,8+1e-9), (5+1*rand(),8), (4,2), (9.5-1e-9,5+1e-9), (1,1.0)])
#test_scan_polygon([(1,0), (1,0), (1,0)])

# ╔═╡ 45c84785-e659-4454-9011-4b3281107498
md"## Find coordinates inside 3D polyhedron"

# ╔═╡ 973a989f-8090-4261-af3b-929f4d418633
begin
	function scan_polyhedron(
	    nodes::Array{Tuple{Float64,Float64,Float64},1},
	    edges::Array{Tuple{Int,Int},1},
	    ε = 1e-9)
	
	    coords = Tuple{Int, Int, Int}[]
	    
	    # find the order of nodes when they are sorted by z, then y, then x
	    nids = sortperm(nodes, by = n -> (n[3], n[2], n[1]))
	    nmin, nmax = nodes[nids[[1,end]]]
	    
	    # reorient the edges such that their first index is the lowest
	    # (according to sort order of nodes)
	    edges = [findfirst(isequal(e[1]), nids) > findfirst(isequal(e[2]), nids) ? (e[2], e[1]) : e for e=edges]
	    
	    # build list of edges starting at the lowest node
	    current = filter(e -> e[1] == nids[1], edges)
	    
	    function distance_to_axis(p) # axis: a+tn
	        a = [nmin...]
	        n = [nmax...] - a
	        n /= (nn=norm(n); nn<ε ? 1 : nn) # do not normalize null vector
	        d = norm((a-p) - dot(a-p, n) * n)
	    end
	    
	    function intersect(e, z)
	        x1, y1, z1 = nodes[e[1]]
	        x2, y2, z2 = nodes[e[2]]
	
	        # return point closer to axis if the points are basically on the same level
	        if abs(z1-z2) < ε
	            d1 = distance_to_axis([x1, y1, z1])
	            d2 = distance_to_axis([x2, y2, z2])
	            return (d1 > d2) ? (x1, y1) : (x2, y2)
	        end
	
	        r = (z-z1) / (z2-z1)
	        (x1 * (1-r) + x2 * r, y1 * (1-r) + y2 * r)
	    end
	    
	    zmin = ceil(Int, nodes[nids[1]][3]-ε)
	    zmax = floor(Int, nodes[nids[end]][3]+ε)
	    
	    # initialize with eager switching
	    i = 2 # index of next node to be passed
	    while (nodes[nids[i]][3]-ε < zmin && i < length(nids))
	        current = [filter(e -> e[2] != nids[i], current);
	                   filter(e -> e[1] == nids[i], edges)]
	        i += 1
	    end
	    
	    for z in zmin:zmax
	        
	        # if a node has been passed, remove edges stopping at i
	        # and add edges starting at i (use lazy switching)
	        while (nodes[nids[i]][3]+ε < z && i < length(nids))
	            current = [filter(e -> e[2] != nids[i], current);
	                       filter(e -> e[1] == nids[i], edges)]
	            i += 1
	        end
	        
	        # intersect edge with z-level to find points of polygon
	        for (x,y) in scan_polygon([intersect(e, z) for e in current])
	            push!(coords, (x, y, z))
	        end
	        
	    end
	    
	    coords
	end
	
	scan_polyhedron(nodes, edges) = 
	    scan_polyhedron([convert(Tuple{Float64,Float64,Float64}, n) for n=nodes], edges)
end

# ╔═╡ d1d6e530-39fa-4600-815b-975095392a86
function test_scan_polyhedron()
    
    f = 2-1e-12
    h = 2*f
    #nodes = [(0,0,f), (0,f,0), (f,0,0),
    #         (h,h,f+h), (h,f+h,h), (f+h,h,h)]
    #nodes = [(0,0,-f), (0,f,0), (f,0,0),
    #         (h,h,h-f), (h,h+f,h), (h+f,h,h)]
    #nodes = [(0,0,0), (f,0,0), (0,f,0),
    #         (0,0,h), (f,0,h), (0,f,h)]
    #nodes = [(0,0,0), (0,f,0), (f,0,0),
    #         (0,0,h), (0,f,h), (f,0,h)]
    nodes = [(2,1,0), (4,1,0), (3,1,1), (2,3,0), (4,3,0), (3,3,1)]
    #nodes = [(2,1,2), (4,1,2), (3,1,1), (2,3,2), (4,3,2), (3,3,1)]
    elems = [(1,2), (2,3), (3,1), (4,5), (5,6), (6,4), (1,4), (2,5), (3,6)]
	fig = Mke.Figure()
	ax = Mke.Axis3(fig[1,1])
    
    for e in elems
        n1 = nodes[e[1]]
        n2 = nodes[e[2]]
        Mke.lines!([n1[1], n2[1]], [n1[2], n2[2]], [n1[3], n2[3]], color=:black)
    end
    
    coords = scan_polyhedron(nodes, elems)
    Mke.scatter!([c[1] for c=coords], [c[2] for c=coords], [c[3] for c=coords], 
		color=:black, markersize=10)

	fig
end

# ╔═╡ 277029cf-135e-465e-9401-f54734caf16e
test_scan_polyhedron()

# ╔═╡ 7d1ca20d-5e04-40ff-a398-b4daba6863ba
md"## Find indices inside band around 3D surface"

# ╔═╡ a9b7c9a1-407f-4ca5-8957-297231aea63c
function generate_cube()
	m = loadmesh("single-cube.stl")
	coords = map(x -> Tuple(Meshes.coordinates(x)), Meshes.vertices(m))
	elems = map(x -> Meshes.indices(x) |> reverse, Meshes.topology(m) |> Meshes.elements)
	(nodes = coords, elems = elems)
end

# ╔═╡ cf3f64fd-249f-44b2-9693-3cb8190f9d46
function test_distance_function(; i3 = nothing)
    
    nx, ny, nz = 300, 300, 20
    lx, ly, lz =  3,  3,  2

    #c = generate_cubes((lx, ly), (1,1), (1,1,1));
    c = generate_cube()
    
    
    # grid for distance function
    x = range(0,lx,nx+1)[1:end-1]
    y = range(0,ly,ny+1)[1:end-1]
    z = range(0,lz,nz+1)
    dx = x[2]-x[1]
    dy = y[2]-y[1]
    dz = z[2]-z[1]
    
    # initialize to infinite distance
    dmax = 0.5
    d = Float64[dmax*2 for x=x, y=y, z=z];
    
    
    function set_coord(coord, val)
        if (1 <= coord[1] <= nx && 1 <= coord[2] <= ny && 1 <= coord[3] <= nz+1)
            if abs(d[coord...]) > abs(val)
                d[coord...] = val
            end
        end
    end
    
    
    for e in c.elems
        
        #println(e)
        n1, n2, n3 = ([c.nodes[ni]...] for ni=e)
        normal = LinearAlgebra.cross(n2-n1, n3-n1)
        
        if abs(normal[2]) < 1e-6 || n1[3] > 0.6 || n2[3] > 0.6 || n3[3] > 0.6
            #continue
        end
        
        normal = normal / LinearAlgebra.norm(normal)
        #println("triangle: ", [n1, n2, n3])
        #println("normal: ", normal)
        
        # add prism shape to distance field
        nd = normal * dmax
        n_prism = [n1-nd, n2-nd, n3-nd, n1+nd, n2+nd, n3+nd]
        e_prism = [(1,4), (2,5), (3,6), (1,2), (2,3), (3,1), (4,5), (5,6), (6,4)]
        
        # transform nodes to index space
        index(pt) = (1 + pt[1]/dx, 1 + pt[2]/dy, 1 + pt[3]/dz)
        value(cd) = [(cd[1]-1)*dx, (cd[2]-1)*dy, (cd[3]-1)*dz]
        i_prism = map(index, n_prism)
        
        #println("n_prism: ", n_prism)
        #println("i_prism: ", i_prism)
        #println("e_prism: ", e_prism)
        
        coords = scan_polyhedron(i_prism, e_prism)
        
        #println(" → found coords: ", coords)
        
        # set coords to one
        for coord in coords
            #xyz = value(coord)
            set_coord(coord, LinearAlgebra.dot(value(coord) - n1, normal))
        end
    end
    
    # plot slice
    xplt=range(x[1]-dx/2, x[end]+dx/2, nx+1)
    yplt=range(y[1]-dy/2, y[end]+dy/2, ny+1)
    #title("z= $(z[izplt])")
    #pcolormesh(xplt, yplt, d[:,:,izplt]')
    #xlim(xplt[[1,end]])
    #ylim(yplt[[1,end]])
    #grid(true)
    #xlabel("x")
    #ylabel("y")
    #colorbar()
	fig = Mke.Figure()
	ax = Mke.Axis(fig[1,1], title="z= $(z[i3])", limits = (xplt[[1,end]]..., yplt[[1,end]]...), xlabel="x", ylabel="y")
	hm = Mke.heatmap!(ax, xplt, yplt, d[:,:,i3])
	Mke.Colorbar(fig[1, 2], hm)

	fig
end

# ╔═╡ 541d710e-8c4d-49fb-be19-7c0d3e9a4c14
test_distance_function(i3 = 8)

# ╔═╡ 9ffeefb3-77ea-42ab-844f-f645b82c87e4
md"## Unused code"

# ╔═╡ e6856840-acce-48bc-a0da-4c3c0d26161b
function scan(pts)
    
    ε = +1e-6
    
    figure(figsize=(6,6))
    ax = gca()
    xlim(-0.25,10.25)
    ylim(-0.25,10.25)
    
    #data = [(i%2==0 && j%2==0 ? 1 : 0) for i=0:20, j=0:20] # checkerboard pattern
    data = [0 for i=0:20, j=0:20]
    xplt = linspace(-0.25,10.25,22)
    yplt = linspace(-0.25,10.25,22)
    
    get(array, i) = array[(n=length(array); i < 1 ? i+n : i>n ? i-n : i)]
    
    x = (xplt[1:end-1]+xplt[2:end])/2
    y = (yplt[1:end-1]+yplt[2:end])/2
    
    ymin, imin = findmin(pt[2] for pt=pts)
    ymax, imax = findmax(pt[2] for pt=pts)
    #ax[:plot]([pts[imin][1]], [pts[imin][2]], "*y")
    #ax[:plot]([pts[imax][1]], [pts[imax][2]], "*g")
    
    wrap(i, n) = i<1 ? i+1 : i>n ? i-n : i
    wrapped_range(i1,i2,n) = i1 <= i2 ? collect(i1:i2) : [collect(i1:n); collect(1:i2)]
    
    n = length(pts)
    if pts[wrap(imin-1,n)][1] < pts[wrap(imin+1,n)][1]
        println("case1")
        lower = reverse(wrapped_range(imax,imin,n))
        upper = wrapped_range(imin,imax,n)
    elseif pts[wrap(imin-1,n)][1] > pts[wrap(imin+1,n)][1]
        println("case2")
        lower = wrapped_range(imin,imax,n)
        upper = reverse(wrapped_range(imax,imin,n))
    end
    
    il = 2
    iu = 2
    println("imin: ", imin, ", ", "imax: ", imax)
    println("lower:", lower)
    println("upper:", upper)
    
    #ax[:plot]([get(pts,i)[1] for i=neighbors], [get(pts,i)[2] for i=neighbors], "--g")
    
    
    
    for iy = find(ymin - ε .< y .< ymax + ε)
        
        yi = y[iy]
        
        while (pts[lower[il]][2] < yi && il < length(lower))
            il += 1
        end
        
        while (pts[upper[iu]][2] < yi && il < length(upper))
            iu += 1
        end
        
        
        
        xl1, yl1 = pts[lower[il-1]]
        xl2, yl2 = pts[lower[il]]
        xu1, yu1 = pts[upper[iu-1]]
        xu2, yu2 = pts[upper[iu]]
        
        xmin = xl1 + (yi-yl1)/(yl2-yl1)*(xl2-xl1)
        xmax = xu1 + (yi-yu1)/(yu2-yu1)*(xu2-xu1)
        println("y = ", yi, ", x ∈ ", [xmin, xmax])
        
        for ix = find(xmin - ε .< x .< xmax + ε)
            data[ix,iy] = 1
        end
        
    end
    
    ax[:pcolormesh](xplt, yplt, data')
    
    ax[:plot]([x for x=x, y=y][:], [y for x=x, y=y][:], ".y")
    ax[:plot]([pt[1] for pt=[pts;pts[1]]], [pt[2] for pt=[pts;pts[1]]], "-r")
    ax[:plot]([pt[1] for pt=pts], [pt[2] for pt=pts], ".k")
    
end

scan(reverse([(1,8), (5+2*rand(),7), (9.5,3), (4,2), (-2.1,3)]))

# ╔═╡ e6029333-77b6-4748-af71-8ddfaff62985
# pts is a selection of points (in order) defining a convex polygon
# xind/yind is a function that takes an x/y-value and returns the corresponding index on a continuous scale
# xval/yval give the x/y value of an index
function scan_polygon2(pts, xind::Function, yind::Function, xval::Function, yval::Function, ε = 1e-6)
    
    # list of indices that are contained in the polygon
    indices = Tuple{Int,Int}[]
    
    # helper function for array of points
    n = length(pts)
    wrap(i) = i<1 ? i+n : i>n ? i-n : i
    function interp_x(p1, p2, y, xrange)
        if abs(p2[2]-p1[2]) < ε
            println("degenerate case!", " p1=",p1, " p2=", p2)
            return (min(xrange[1], p1[1], p2[1]), max(xrange[2], p1[1], p2[1]))
        end
        x = p1[1] + (p2[1]-p1[1]) * (y-p1[2])/(p2[2]-p1[2])
        (min(xrange[1], x), max(xrange[2], x))
    end
    
    # find y-extent of polygon
    ymin, imin = findmin(pt[2] for pt=pts)
    ymax, imax = findmax(pt[2] for pt=pts)
    
    # indices for going left & right along the polygon
    il = wrap(imin-1)
    ir = wrap(imin+1)
    
    for iy = ceil(Int, yind(ymin-ε)):floor(Int, yind(ymax+ε))
        
        y = yval(iy)
        
        # advance currently active line segments, not going past last point
        # TODO: add epsilon for switch at correct moment
        while (pts[il][2] < y && il != imax)
            il = wrap(il-1)
            println("advancing il to $(il) at y=$(y)")
        end
        while (pts[ir][2] < y && ir != imax)
            ir = wrap(ir+1)
            println("advancing ir to $(ir) at y=$(y)")
        end
        
        # find x-values for current y-value on line segments
        xrange = [Inf, -Inf]
        xrange = interp_x(pts[il], pts[wrap(il+1)], y, xrange)
        xrange = interp_x(pts[wrap(ir-1)], pts[ir], y, xrange)
        
        println("y = ", y, ", x ∈ ", xrange)
        
        for ix = ceil(Int, xind(xrange[1]-ε)):floor(Int, xind(xrange[2]+ε))
            push!(indices, (ix,iy))
        end
        
    end
    
    return indices
end

# prepare input
data = [0 for i=0:20, j=0:20]
pts = ([(1,8), (5+2*rand(),8), (9.6,2), (4,2), (1,3.0)])
xval(i) = (i-1)/2
yval(j) = (j-1)/2
xind(x) = 1 + 2*x
yind(y) = 1 + 2*y

# compute output & add to data
indices = scan_polygon2(pts, xind, yind, xval, yval)

for i in indices
    #println("i=", i, " of type", typeof(i))
    data[i...] = 1
end

# prepare figure
figure(figsize=(6,6))
ax = gca()
xlim(-0.25,10.25)
ylim(-0.25,10.25)

# draw pixels
xplt = linspace(-0.25,10.25,22)
yplt = linspace(-0.25,10.25,22)
xvals = (xplt[1:end-1]+xplt[2:end])/2
yvals = (yplt[1:end-1]+yplt[2:end])/2
ax[:pcolormesh](xplt, yplt, data')
ax[:plot]([x for x=xvals, y=yvals][:], [y for x=xvals, y=yvals][:], ".w", ms=1)

# draw polygon
ax[:plot]([pt[1] for pt=[pts;pts[1]]], [pt[2] for pt=[pts;pts[1]]], "-r")
ax[:plot]([pt[1] for pt=pts], [pt[2] for pt=pts], ".k")

pygui(false)

# ╔═╡ cf1669a2-1096-4015-814e-7b423b8b051e
function intersect_z(p1, p2, z)
    ε = +1e-6
    
    # return first value if the points are basically on the same level
    if abs(p1[3]-p2[3]) < ε
        return (p1[1], p1[2])
    end
    
    r = (z-p1[3]) / (p2[3]-p1[3])
    (p1[1] * (1-r) + p2[1] * r,
     p1[2] * (1-r) + p2[2] * r)
end

function intersect_prism(n1, n2, n3, v, xvals, yvals, zvals)
    assert(n1[3] <= n2[3] <= n3[3])
    assert(v[3] >= 0)
    ε = +1e-6
    
    N1 = [n1[i]+v[i] for i=1:3]
    N2 = [n2[i]+v[i] for i=1:3]
    N3 = [n3[i]+v[i] for i=1:3]
    
    println("n1,n2,n3=", n1, n2, n3)
    println("N1,N2,N3=", N1, N2, N3)
    println("v=", v)
    
    layer_edges = [(n1,n2), (n1,n3), (n1,N1)]
    
    fig = figure()
    ax = gca(projection="3d")
    xlabel("x")
    ylabel("y")
    
    # prepare functions for switching between values & indices for x & y
    x0, y0 = xvals[1], yvals[1]
    dx, dy = xvals[2]-xvals[1], yvals[2]-yvals[1]
    xval(i) = x0 + (i-1) * dx
    yval(j) = y0 + (j-1) * dy
    xind(x) = 1 + (x-x0) / dx
    yind(y) = 1 + (y-y0) / dy
    
    # plot dot grid
    #ax[:plot]([x for x=xvals, y=yvals, z=zvals][:], [y for x=xvals, y=yvals, z=zvals][:],
    #    [z for x=xvals, y=yvals, z=zvals][:], ".k", ms=1)
    
    
    for z in zvals[n1[3] - ε .< zvals .< N3[3] + ε]
        
        while true
            replace = [e[2][3] < z + 2ε for e = layer_edges]
            println(replace)
            if any(replace)
                nreplace = find(replace)[1]
                replacenode = layer_edges[nreplace][2]
                replacement = if replacenode==n2
                    println("passed n2 node at z=", z)
                    [(n2,N2), (n2,n3)] # order matters!
                elseif replacenode==n3
                    println("passed n3 node at z=", z)
                    [(n3,N3)]
                elseif replacenode==N1
                    println("passed N1 node at z=", z)
                    [(N1,N3),(N1,N2)] # order matters!
                elseif replacenode==N2
                    println("passed N2 node at z=", z)
                    [(N2,N3)]
                else # N3 has been surpassed
                    println("done with prism at z=", z)
                    return
                end
                layer_edges = filter(e->e[2]!=replacenode,
                    [layer_edges[1:nreplace-1]; replacement...; layer_edges[nreplace+1:end]])
                println("layer_edges: ", layer_edges)
                continue
            end
            break
        end
            
        
        pts = [intersect_z(e..., z) for e in layer_edges]
        i_pgn = scan_polygon(pts, xind, yind, xval, yval, ε)
        println("scan polygon for layer z=",z, " with polygon ", pts)
        
        
        # plot polygon & points
        ax[:plot]([pt[1] for pt=[pts;pts[1]]], [pt[2] for pt=[pts;pts[1]]], [z for pt=[pts;pts[1]]], "-r")
        ax[:plot]([xval(i[1]) for i=i_pgn], [yval(i[2]) for i=i_pgn], [z for i=i_pgn], ".g", ms=4)
    end
end

pygui(true)

#=
xv = linspace(-15, 5, 41)
yv = linspace(-15, 5, 41)
zv = linspace(0,10,41)

# case n1 < n2 < n3 < N1 < N2 < N3
v01 = rand(3)
v12=[1,0,2]
v13=[0,1,3]
v1n = cross(v12, v13) * 5
assert(v1n[3] > 3)
intersect_prism(v01, v01+v12, v01+v13, v1n, xv, yv, zv)
    =#

xv = linspace(0, 1.5, 21)
yv = linspace(0, 1.5, 21)
zv = linspace(0,20,41)
    v01 = [0,0,3]
    v12 = [1,0,0]
    v13 = [0,1,0]
    v1n = [0,0,9]
intersect_prism(v01, v01+v12, v01+v13, v1n, xv, yv, zv)

# ╔═╡ b4f5db9d-c401-4423-b419-a0e82c9bcfe6
# - let’s say, we always go from the lower side to the higher side
# - at most, the prism can be horizontal (limiting case)
# - the normal is either positive, or horizontal
# - note: the triangle can still be oriented in any direction, no restriction on surface
# - we start with the lowest point of the triangle
# - since the normal is upwards, the highest point is on the other triangle
# - we can order the remaining four points by z-index
# - therefore, we can have zero, one, two, three, or four inner points below the z-index
# - all three points get an equal increase of z, so the lowest point on the first triangle
#   is the same as the lowest on the second one
# - zero points: polygon is a triangle, all three vectors from lowest point intersect with plane
# - one point: could be from triangle or from following the normal
#   - triangle: 

# ╔═╡ fd662a4f-768e-4f24-b4d4-3b1cf0ed3866
z0, dz = z[1], z[2]-z[1]
ε = 1e-6
dmax = 0.5

function add_polygon(n1, n2, n3, normal, iz)
    zlvl = z[iz]
    o  = [n1...] # origin of prism (0->n1)
    v1 = [n2...] - [n1...] # vector (n1->n2)
    v2 = [n3...] - [n2...] # vector n2 -> n3
    v3 = normal # only one side for now
    
end

for e in c.elems
    n1, n2, n3 = (c.nodes[ni] for ni=e)
    normal = cross([n2[i]-n1[i] for i=1:3], [n3[i]-n2[i] for i=1:3])
    normal = normal / norm(normal) * dmax
    #if (abs(normal[3]) > ε) # ignore vertical triangles for now
    zmin = min(n1[3], n2[3], n3[3]) + min(normal[3],0) - ε
    zmax = max(n1[3], n2[3], n3[3]) + max(normal[3],0) + ε
    izmin = 1 + max(0, ceil(Int, (zmin-z0)/dz))
    izmax = 1 + floor(Int, (zmax-z0)/dz)
    add_polygon(n1, n2, n3, normal, izmin)

    println("norm: ", normal, ", zrange: ", (zmin, zmax), ", izrange: ", (izmin, izmax))
end

# ╔═╡ Cell order:
# ╟─9929f2be-277a-11ed-375a-f18a1b52c0ca
# ╟─87aca29f-29df-44ba-9630-cc97e927a339
# ╟─3671d217-c327-431f-9bfd-0f3e7fd4a01a
# ╟─fc1826d8-f6ad-443f-bf65-2d5def8974a7
# ╟─c4697b84-9cbd-4178-9aaf-72527bddde67
# ╠═f33325a0-ca9e-4366-bc37-91dc414cad25
# ╠═34b3f98c-1053-4aa7-ae0f-4aff49998c82
# ╟─82d773d2-a328-41ca-b5c6-1836a393d6da
# ╟─07d89b37-ffc9-49c1-9f7c-59e26f256fcd
# ╟─a45d21c3-ea3b-4d11-8a4a-be519f1d7719
# ╟─17854bbc-543c-42bd-a126-943edda8c75a
# ╟─2d4b6ea3-decb-4f6a-939d-1eafe956828e
# ╟─3772c9e9-2db1-4978-a011-de19799f4e98
# ╟─cb2c81fe-1987-42c9-aec8-badbbd0b1a0a
# ╟─0d7b073f-90a3-4588-b4d1-5cc0b14e2bc5
# ╟─409498df-51ed-49e8-85a9-efcdee401255
# ╠═ef3d6f6e-327f-45b8-92c6-423dd52bac49
# ╠═e0415fec-1d32-458d-ae55-87251486cbcb
# ╠═01843255-f9a3-49a7-a989-b15b2e6932de
# ╟─3c555464-59e2-4c98-bdfe-f64891f6b569
# ╟─bf425675-9ee5-43de-a67a-b9098a751187
# ╟─342ed4c3-ffcc-4a75-887d-5c9c6f8318d6
# ╠═73dfeea2-3faf-458c-ad41-7d592a581479
# ╠═1f405903-9e56-408c-aea6-1f1328bc58c9
# ╠═365fd982-261a-4c11-85d3-0f54b6d4865a
# ╠═3d03fcd5-a804-452e-b409-8d7317df014f
# ╠═1b06df66-f0db-40c6-a8ad-79e855c883e0
# ╠═3bd42cac-c944-44c3-a8d9-e4eb07d059d3
# ╠═a1588120-35c9-4719-8bca-757f56cc288a
# ╟─ff1d3625-60e0-41dd-b889-0d2606378568
# ╟─a07eeeed-1d9d-4574-bad3-fabc4f0b0af9
# ╟─d4065db9-6b48-4c1e-ae33-ed793fc91c28
# ╠═762e6049-bf06-437c-b990-590a2d025021
# ╠═f809754a-ea9f-4062-8a94-85e64b42c3a2
# ╠═f415eab0-b090-495c-a79c-64a8ad4f3232
# ╟─45c84785-e659-4454-9011-4b3281107498
# ╠═277029cf-135e-465e-9401-f54734caf16e
# ╟─973a989f-8090-4261-af3b-929f4d418633
# ╠═d1d6e530-39fa-4600-815b-975095392a86
# ╟─7d1ca20d-5e04-40ff-a398-b4daba6863ba
# ╠═541d710e-8c4d-49fb-be19-7c0d3e9a4c14
# ╠═a9b7c9a1-407f-4ca5-8957-297231aea63c
# ╠═cf3f64fd-249f-44b2-9693-3cb8190f9d46
# ╟─9ffeefb3-77ea-42ab-844f-f645b82c87e4
# ╠═e6856840-acce-48bc-a0da-4c3c0d26161b
# ╠═e6029333-77b6-4748-af71-8ddfaff62985
# ╠═cf1669a2-1096-4015-814e-7b423b8b051e
# ╠═b4f5db9d-c401-4423-b419-a0e82c9bcfe6
# ╠═fd662a4f-768e-4f24-b4d4-3b1cf0ed3866
