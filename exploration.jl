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

# ╔═╡ a96a6610-2199-414b-a97b-2f3c4a8d33d2
function loadmesh(path; flip = false)
	mesh = FileIO.load(path)
	points = [Tuple(p) for p in Set(mesh.position)]
	indices = Dict(p => i for (i, p) in enumerate(points))
	connectivities = map(mesh) do el
		(i,j,k) = (indices[Tuple(p)] for p in el)
		Meshes.connect(flip ? (i,k,j) : (i,j,k))
	end
    Meshes.SimpleMesh(Meshes.Point3.(points), connectivities)
end

# ╔═╡ f33325a0-ca9e-4366-bc37-91dc414cad25
let mesh = loadmesh("single-cube.stl")
	MeshViz.viz(mesh, showfacets=true, color=1:Meshes.nelements(mesh))
end

# ╔═╡ 4ca47d96-e749-4b06-9b1e-8fa2222a0936
md"## Rasterization Algorithms"

# ╔═╡ 908217b7-39f8-42bd-870b-92a2e872b7a7
begin
	abstract type Rasterization end
	struct BoundingBox <: Rasterization end
	struct EdgeFunctions{T} <: Rasterization
		tol::T
		EdgeFunctions(tol = 1e-9) = new{typeof(tol)}(tol)
	end
end

# ╔═╡ 5a1e55bb-7a96-46c4-85d0-cf40606860bf
md"### Extension of Polyhedra"

# ╔═╡ 73dfeea2-3faf-458c-ad41-7d592a581479
struct Prism{Dim,T,V<:AbstractVector{Meshes.Point{Dim,T}}} <: Meshes.Polyhedron{Dim,T}
	base::V # should be oriented such that the normal and the side are colinear
    side::Meshes.Vec{Dim,T}

	function Prism(base::V, side, fix=true) where {Dim,T,V<:AbstractVector{Meshes.Point{Dim,T}}}
		if fix && sign(Meshes.dot(Meshes.cross(base[2]-base[1], base[3]-base[1]), side)) != 1
			base = reverse(base)
		end
		new{Dim,T,V}(base, side)
	end
end

# ╔═╡ c3ce016a-b6bb-4bf6-aa48-4f5da7818c77
begin
	connections(p::Prism, rank) = if rank == 0
		((1,), (2,), (3,), (4,), (5,), (6,))
	elseif rank == 1
		((1,2), (2,3), (3,1), (1,4), (2,5), (3,6), (4,5), (5,6), (6,4))
	elseif rank == 2
		((3,2,1),(5,4,1,2),(3,6,5,2),(1,4,6,3),(5,6,4))
	elseif rank == 3
		((1, 2, 3, 4, 5, 6))
	else
		error("Invalid rank '$rank' for prism")
	end

	# NOTE: creating a boundary allocates a vector for the connections
	Meshes.boundary(p::Prism) = Meshes.SimpleMesh(Meshes.vertices(p), collect(Meshes.connect.(connections(p, 2))))

	Meshes.nvertices(p::Prism) = 6
	Meshes.vertices(p::Prism) = Meshes.Vec(p.base[1], p.base[2], p.base[3],
		p.base[1] + p.side, p.base[2] + p.side, p.base[3] + p.side)

	nedges(p::Meshes.Polytope) = Meshes.nfaces(Meshes.boundary(p), 1)
	nedges(p::Prism) = 9
	edges(p::Meshes.Polytope) = Meshes.faces(Meshes.boundary(p), 1)
	edges(p::Prism) = let v = Meshes.vertices(p), c = connections(p, 1)
		(Meshes.materialize(Meshes.connect(c), v) for c in c)
	end

	Meshes.nfacets(p::Meshes.Polytope) = Meshes.nfaces(Meshes.boundary(p), Meshes.paramdim(p)-1)
	Meshes.nfacets(p::Prism) = 5
	Meshes.facets(p::Meshes.Polytope) = Meshes.faces(Meshes.boundary(p), Meshes.paramdim(p)-1)
	Meshes.facets(p::Prism) = let v = Meshes.vertices(p), c = connections(p, 2)
		(Meshes.materialize(Meshes.connect(c), v) for c in c)
	end

	Meshes.nfaces(p::Meshes.Polytope, rank) = if rank == 0
		Meshes.nvertices(p)
	elseif rank == 1
		nedges(p)
	elseif rank == Meshes.paramdim(p) - 1
		Meshes.nfacets(p)
	else
		Meshes.nfaces(Meshes.boundary(p), rank)
	end
	Meshes.faces(p::Meshes.Polytope, rank) = if rank == 0
		Meshes.vertices(p)
	elseif rank == 1
		edges(p)
	elseif rank == Meshes.paramdim(p) - 1
		Meshes.facets(p)
	else
		Meshes.faces(Meshes.boundary(p), 1)
	end
end

# ╔═╡ 34b3f98c-1053-4aa7-ae0f-4aff49998c82
let dims = (3, 3, 3), origin = (0.25, 0.25, 0.25), spacing = (0.5, 0.5, 0.5)
	gd = Meshes.CartesianGrid(dims, origin, spacing)
	gd = Meshes.vertices(gd)
	MeshViz.viz(gd, size=10)
end

# ╔═╡ 8b9a30c3-f4ab-4367-b2a0-a398da18cd1b
begin

	struct Hyperplane{Dim,T} <: Meshes.Primitive{Dim,T}
		p::Meshes.Point{Dim,T}
		n::Meshes.Vec{Dim,T}
	end

	"""
		hyperplanes(polytope)

	Return the hyperplanes in which the boundaries of a polytope live, oriented such that the plane normals point out of the polytope.
	"""
	function hyperplanes(pt::Meshes.Polytope{Dim,Dim,T}) where {Dim,T}
		error("Not implemented")
	end

	function hyperplanes(pg::Meshes.Ngon{N,2,T}) where {N,T}
		v = Meshes.vertices(pg)
		ntuple(N) do i
			p = v[i]
			dp = p - (i == 1 ? v[end] : v[i-1])
			n = Meshes.Vec(dp[2], -dp[1])
			Hyperplane(p, n)
		end
	end

	function hyperplanes(ph::Meshes.Polyhedron{3,T}) where {T}
		# TODO: this only works correctly for some polyhedra
		Iterators.map(Meshes.facets(ph)) do facet
			xs = Meshes.vertices(facet)
			p::Meshes.Point{3,T} = xs[1]
			dp1::Meshes.Vec{3,T} = xs[2] - p
			dp2::Meshes.Vec{3,T} = xs[3] - p
			Hyperplane(p, Meshes.cross(dp1, dp2))
		end
	end

	function hyperplanes(pr::Prism{3,T}) where {T}
		v1 = pr.base[2] - pr.base[1]
		v2 = pr.base[3] - pr.base[2]
		v3 = pr.base[1] - pr.base[3]
		v4 = pr.side
		vs = (v1, v2, v3, v4)
		ntuple(5) do i
			p = isodd(i) ? pr.base[1] : pr.base[2] + pr.side
			dp1 = vs[(1,2,3,1,1)[i]]
			dp2 = vs[(4,4,4,2,3)[i]]
			Hyperplane(p, Meshes.cross(dp1, dp2))
		end
	end
end

# ╔═╡ a45d21c3-ea3b-4d11-8a4a-be519f1d7719
begin

	function rasterize(box::Meshes.Box{Dim,T},
					   grid::Meshes.CartesianGrid{Dim,T}) where {Dim,T}
		dx = Meshes.spacing(grid)
		imin = (minimum(box) - minimum(grid)) ./ dx .+ 1 |>
			i -> floor.(Int, i) |>
			i -> max.(i, 1)
		imax = (maximum(box) - minimum(grid)) ./ dx |>
			i -> ceil.(Int, i) |>
			i -> min.(i, size(grid))
		inds = ntuple(i -> imin[i]:imax[i], Dim)
		CartesianIndices(inds)
	end

	rasterize(geom, grid) = rasterize(geom, grid, BoundingBox())

	# test if point in polygon for all points in bounding box
	function rasterize(geom::Meshes.Geometry{Dim,T},
					   grid::Meshes.CartesianGrid{Dim,T},
		               ::BoundingBox) where {Dim,T}

		Δgd = Meshes.spacing(grid) |> Meshes.Vec
		xref = minimum(grid) - 0.5 * Δgd

		inds = rasterize(Meshes.boundingbox(geom), grid)
		Iterators.filter(inds) do i
			pt = xref + Δgd .* i.I
			pt ∈ geom
		end
	end

	function rasterize(polytope::Meshes.Polytope{Dim,Dim,T},
					   grid::Meshes.CartesianGrid{Dim,T},
					   method::EdgeFunctions) where {Dim,T}
		E = edgefunctions(polytope, grid)
		tol = convert(T, method.tol)
		inds = rasterize(Meshes.boundingbox(polytope), grid)
		Iterators.filter(inds) do ind
			all(apply_edgefunction(Ei, ind) < tol for Ei in E)
		end
	end

	function edgefunctions(polytope::Meshes.Polytope{Dim,Dim,T},
				grid::Meshes.CartesianGrid{Dim,T},
				pos = T(0.5)) where {Dim,T}
		dx = Meshes.spacing(grid) |> Meshes.Vec
		xref = minimum(grid) + (pos .- 1) .* dx
		map(hyperplanes(polytope)) do hp
			# as the facets should have outward-pointing normals, the edge functions
			# as defined here should be negative for points inside the polyhedron
			dE = hp.n .* dx
			Eref = Meshes.dot(xref - hp.p, hp.n)
			Meshes.Vec(dE..., Eref)
		end
	end

	# explicit multiplication appears to run faster than dot product
	@inline apply_edgefunction(Ei::Meshes.Vec{3}, ind::CartesianIndex{2}) =
		Ei[1] * ind[1] + Ei[2] * ind[2] + Ei[3]
	@inline apply_edgefunction(Ei::Meshes.Vec{4}, ind::CartesianIndex{3}) =
		Ei[1] * ind[1] + Ei[2] * ind[2] + Ei[3] * ind[3] + Ei[4]

	rasterize
end

# ╔═╡ 6969c27d-6609-4f30-94fd-5d1c9a4fdc91
# ╠═╡ disabled = true
#=╠═╡
begin
	function rasterize3d(poly, grid, ::Meshes.GrahamScan)

		dx, dy, dz = Meshes.spacing(grid)
		z0 = Meshes.coordinates(minimum(grid))[3] - 0.5 * dz

		dims2d = size(grid)[1:2]
		origin2d = Meshes.Point(Meshes.coordinates(minimum(grid))[1:2])
		spacing2d = Meshes.SVector(dx, dy)
		grid2d = Meshes.CartesianGrid(dims2d, origin2d, spacing2d)

		segs = Meshes.faces(poly, 1)
		inds = rasterize(Meshes.boundingbox(poly), grid)

		Iterators.map(inds.indices[3]) do iz

			# find intersection points with z-plane
			pl = Meshes.Plane((0.,0.,z0 + iz * dz), (1.,0.,0.), (0.,1.,0.))
			pts = segs |>
				sgs -> Iterators.map(x->Meshes.intersect(pl, x), sgs) |>
				pts -> Iterators.filter(x->!isnothing(x), pts) |>
				pts -> Iterators.map(x->Meshes.Point2(Meshes.coordinates(x)[1:2]...), pts) |>
				collect

			# skip layer if no intersections
			isempty(pts) && return ()

			# construct polygon
			ng = Meshes.hull(Meshes.PointSet(pts), Meshes.GrahamScan()) |>
				Meshes.vertices |> Meshes.Ngon

			Iterators.map(rasterize(ng, grid2d, EdgeFunctions())) do ind
				CartesianIndex(ind.I..., iz)
			end
		end |> Iterators.flatten
	end
end
  ╠═╡ =#

# ╔═╡ 82d773d2-a328-41ca-b5c6-1836a393d6da
md"## Rasterize 2D Polygons"

# ╔═╡ 7999f361-6a5b-4712-b7cb-3904fbecd108
let
	N = (2, 3)
	gd = Meshes.CartesianGrid(N, Meshes.Point(0, 0), (10 ./ N))
	dx = Meshes.spacing(gd)
	minimum(gd) + Meshes.Vec(dx)
end

# ╔═╡ 07d89b37-ffc9-49c1-9f7c-59e26f256fcd
let
	#pg = Meshes.Ngon([(1,8+1e-9), (5+1*rand(),8), (4,2), (9.5-1e-9,5+1e-9), (1,1.0)])

	# set up polygon
	pts = Meshes.PointSet([(1,8+1e-9), (5+1*rand(),8), (4,2), (9.5-1e-9,5+1e-9), (1,1.0)])
	pg = Meshes.hull(pts, Meshes.GrahamScan()) |> Meshes.vertices |> collect |> Meshes.Ngon

	N = (32, 32)

	gd = Meshes.CartesianGrid(N, Meshes.Point(0, 0), (10 ./ N))
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
	 broken = ()), # TraverseBoundingBox may or may not pass
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
			@test [x.I for x in rasterize(ng, gd, method())] ==
				collect(case.result) broken=broken
		end
	end
end

# ╔═╡ 2d4b6ea3-decb-4f6a-939d-1eafe956828e
@testset "Rasterization" begin
	runtests(test_cases, (BoundingBox, EdgeFunctions))
end;

# ╔═╡ 0d7b073f-90a3-4588-b4d1-5cc0b14e2bc5
md"### Benchmarks"

# ╔═╡ 409498df-51ed-49e8-85a9-efcdee401255
function runbenchmarks2d(method)
	ngon = Meshes.Ngon([(1.2,1.7), (1.8,1.8), (1.3,2.1), (1.1,1.9)])
	grid = Meshes.CartesianGrid((256,256), (0.,0.), (4π/256, 2π/256))
	@benchmark maximum(rasterize($ngon, $grid, $method))
end

# ╔═╡ e0415fec-1d32-458d-ae55-87251486cbcb
# ╠═╡ disabled = true
#=╠═╡
runbenchmarks2d(BoundingBox())
  ╠═╡ =#

# ╔═╡ 01843255-f9a3-49a7-a989-b15b2e6932de
# ╠═╡ disabled = true
#=╠═╡
runbenchmarks2d(EdgeFunctions())
  ╠═╡ =#

# ╔═╡ 3c555464-59e2-4c98-bdfe-f64891f6b569
md"### Allocations"

# ╔═╡ bf425675-9ee5-43de-a67a-b9098a751187
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

# ╔═╡ 8148d1c8-8fea-4a7d-a1d4-23a256a10b8a
let
	ngon = Meshes.Ngon([(1.2,1.7), (1.8,1.8), (1.3,2.1), (1.1,1.9)])
	grid = Meshes.CartesianGrid((256,256), (0.,0.), (4π/256, 2π/256))
	#method = TraverseBoundingBox()
	method = EdgeFunctions()
	function f(N)
		imax = CartesianIndex(1,1)
		for i=1:N
			imax = maximum(rasterize(ngon, grid, method), init=imax)
		end
		imax
	end
	f(1000)
	Profile.clear()
	Profile.@profile f(100000)
	PProf.pprof(web=false)
end

# ╔═╡ 342ed4c3-ffcc-4a75-887d-5c9c6f8318d6
md"## Rasterize 3D Prisms"

# ╔═╡ 1dd29d36-dfd3-4bea-a29c-b190c5713575
let
    f = 2-1e-12
    h = 2*f
    #pr = Prism(Meshes.Point.([(0,0,f), (0,f,0), (f,0,0)]), (h,h,h))
    #pr = Prism(Meshes.Point.([(0,0,-f), (0,f,0), (f,0,0)]), (h,h,h))
    #pr = Prism(Meshes.Point.([(0,0,0), (0,f,0), (f,0,0)]), (0,0,h))
	pr = Prism(Meshes.Point.([(2,1,0), (4,1,0), (3,1,1)]), (0,2,0))
	#pr = Prism(Meshes.Point.([(2,1,2), (4,1,2), (3,1,1)]), (0,2,0))

	gd = Meshes.CartesianGrid((0.,0.,0.), (5.,5.,5.), dims=(10,10,10))
	#gd = Meshes.CartesianGrid((0,0,0), (4,4,4), dims=(4,4,4))
	#gd = Meshes.CartesianGrid((-0.5,-0.5,-0.5), (4.5,4.5,4.5), dims=(5,5,5))

	#method = Meshes.GrahamScan()
	method = EdgeFunctions()

	pts = map(rasterize(pr, gd, method)) do i
		Meshes.centroid(gd[LinearIndices(size(gd))[i]])
	end

	fig = Mke.Figure(resolution = (400, 400))
	ax = Mke.Axis3(fig[1,1])
	MeshViz.viz!(Meshes.PointSet(Meshes.vertices(pr)), color=:black)
	isempty(pts) || MeshViz.viz!(pts, color=:yellow)
	MeshViz.viz!(pr, alpha=0.25, transparency=true)
	fig
end

# ╔═╡ c49f0f3d-043c-4dc4-aab8-7fa7b5a23dda
md"### Checks & Benchmarks"

# ╔═╡ 70f8f6cb-47d2-4254-8032-f93a713f8e74
let
	prism = Prism(Meshes.Point.([(1.2,1.7,0.1), (1.8,1.8,0.1), (1.3,2.1,0.2)]), (0.1,-0.1,0.5))
	grid = Meshes.CartesianGrid((256,256,256), (0.,0.,0.), (4π/256, 2π/256, 1/256))
	@code_warntype rasterize(prism, grid, EdgeFunctions())
end

# ╔═╡ ad748fba-054c-4d74-b589-24dbc0d8afd1
#=╠═╡
let
	prism = Prism(Meshes.Point.([(1.2,1.7,0.1), (1.8,1.8,0.1), (1.3,2.1,0.2)]), (0.1,-0.1,-0.5))
	grid = Meshes.CartesianGrid((256,256,256), (0.,0.,0.), (4π/256, 2π/256, 1/256))
	rg = rasterize3d(prism, grid, Meshes.GrahamScan()) |> collect
	re = rasterize3d(prism, grid, EdgeFunctions()) |> collect
	@test rg == re
end
  ╠═╡ =#

# ╔═╡ 8ffe0ca3-07bd-46ff-b525-8890acf98b68
function runbenchmarks3d(method)
	prism = Prism(Meshes.Point.([(1.2,1.7,0.1), (1.8,1.8,0.1), (1.3,2.1,0.2)]), (0.1,-0.1,0.5))
	grid = Meshes.CartesianGrid((256,256,256), (0.,0.,0.), (4π/256, 2π/256, 1/256))
	@benchmark maximum(rasterize($prism, $grid, $method))
end

# ╔═╡ a1588120-35c9-4719-8bca-757f56cc288a
# ╠═╡ disabled = true
#=╠═╡
runbenchmarks3d(EdgeFunctions())
  ╠═╡ =#

# ╔═╡ ba3a5975-631f-44d3-9d7f-96c89927b7a2
runbenchmarks3d(EdgeFunctions())

# ╔═╡ 59335790-2400-4e6b-aff5-5f3837503033
# ╠═╡ disabled = true
#=╠═╡
let
	prism = Prism(Meshes.Point.([(1.2,1.7,0.1), (1.8,1.8,0.1), (1.3,2.1,0.2)]), (0.1,-0.1,0.5))
	grid = Meshes.CartesianGrid((256,256,256), (0.,0.,0.), (4π/256, 2π/256, 1/256))
	method = EdgeFunctions()
	function f(N)
		imax = CartesianIndex(1,1,1)
		for i=1:N
			imax = maximum(rasterize(prism, grid, method), init=imax)
		end
		imax
	end
	f(1000)
	Profile.clear()
	Profile.@profile f(1000)
	PProf.pprof(web=false)
end
  ╠═╡ =#

# ╔═╡ 86244339-1362-48b0-8342-d4a6528138e3
# ╠═╡ disabled = true
#=╠═╡
let
	prism = Prism(Meshes.Point.([(1.2,1.7,0.1), (1.8,1.8,0.1), (1.3,2.1,0.2)]), (0.1,-0.1,0.5))
	grid = Meshes.CartesianGrid((256,256,256), (0.,0.,0.), (4π/256, 2π/256, 1/256))
	#method = TraverseBoundingBox()
	method = EdgeFunctions()
	maximum(rasterize(prism, grid, method), init = CartesianIndex(1,1,1))
	Profile.Allocs.clear()
	Profile.Allocs.start(sample_rate=1.0)
	maximum(rasterize(prism, grid, method), init = CartesianIndex(1,1,1))
	Profile.Allocs.stop()
	PProf.Allocs.pprof(web=false)
end
  ╠═╡ =#

# ╔═╡ fccdb7f1-4ca3-47ec-85a9-8d6b0b2e23c8
md"## Rasterize Balls"

# ╔═╡ e24d37db-79ea-45fc-a5f4-2b372b30f5a8
let
	b = Meshes.Ball((2,2,2), 1.5)

	gd = Meshes.CartesianGrid((0.,0.,0.), (5.,5.,5.), dims=(10,10,10))

	pts = map(rasterize(b, gd)) do i
		Meshes.centroid(gd[LinearIndices(size(gd))[i]])
	end

	fig = Mke.Figure(resolution = (400, 400))
	ax = Mke.Axis3(fig[1,1])
	MeshViz.viz!(Meshes.PointSet(Meshes.vertices(gd)), color=:black, size=3)
	isempty(pts) || MeshViz.viz!(pts, color=:yellow)
	MeshViz.viz!(b, alpha=0.25, transparency=true)
	fig
end

# ╔═╡ 35e13f7a-cde4-430c-9fe8-51c96c2e5af4
function runbenchmarks_ball()
	b = Meshes.Ball((2,2,2), 1.5)
	gd = Meshes.CartesianGrid((0.,0.,0.), (5.,5.,5.), dims=(10,10,10))
	@benchmark maximum(rasterize($b, $gd))
end

# ╔═╡ 86e0b953-924e-44d2-acd4-1515e732eb6e
runbenchmarks_ball()

# ╔═╡ 9069a1a1-31b2-4e5a-80fd-a410475f4d88
md"## Distance Field"

# ╔═╡ 3e13d529-2114-4495-a1eb-7669f2a068c7
function plane_extrusions(ds, ns, mesh, grid, dmax; verbose = true)
	# process distances to triangle surface
	lininds = LinearIndices(size(grid))
	for el in mesh
		n = Meshes.normal(el)
		v = Meshes.vertices(el)
		prism = Prism(v .- (dmax * n,), 2 * dmax * n)
		for i in rasterize(prism, grid, EdgeFunctions())
			x = Meshes.centroid(grid[lininds[i]])
			d = Meshes.dot(x - v[1], n)
			d_old = ds[i]
			if ismissing(d_old) || abs(d) < abs(d_old)
				ds[i] = d
			end
		end
	end
end

# ╔═╡ feaefa19-fb04-42c5-bdd7-b095542b9cd2
function edge_extrusions(ds, ns, mesh, grid, dmax; verbose = true)

	lininds = LinearIndices(size(grid))

	# process distances to triangle edges
	topo = convert(Meshes.HalfEdgeTopology, Meshes.topology(mesh))
	e2v = Meshes.Boundary{1,0}(topo)
	e2t = Meshes.Coboundary{1,2}(topo)
	flat_edges = 0
	border_edges = 0
	for ie in 1:Meshes.nfacets(topo)

		# determine connections of edges
		neighbors = e2t(ie)
		if length(neighbors) == 1
			border_edges += 1
			continue # border edges are skipped
		end
		nb1, nb2 = neighbors
		n1, n2 = Meshes.normal.((mesh[nb1], mesh[nb2]))
		p1, p2 = Meshes.vertices(mesh)[e2v(ie)]

		# determine if edge is convex or concave
		o = sign(Meshes.dot(Meshes.cross(p2-p1, n1), n2)) # positive for convex edge
		if o == 0
			flat_edges += 1
			continue # skip flat edges
		end

		# mean normal at edge
		nmid = n1 + n2
		lmid = Meshes.norm(nmid)
		lmid > 1e-3 || error("Degenerate mesh: acute edge angle")
		nmid /= lmid

		# length of prism sides required that distance is at least dmax everywhere
		lside = dmax / Meshes.dot(nmid, n1)

		# define prism for edge extrusion
		base = Meshes.Vec(p1, p1 + lside * o * n1, p1 + lside * o * n2)
		side = p2 - p1
		prism = Prism(base, side)

		# compute distance for points inside prism
		for i in rasterize(prism, grid, EdgeFunctions())
			x = Meshes.centroid(grid[lininds[i]])
			d = o * Meshes.norm(Meshes.cross(x - p1, side)) / Meshes.norm(side)
			d_old = ds[i]
			if ismissing(d_old) || abs(d) < abs(d_old)
				if abs(d) <= dmax
					ds[i] = d
				end
			end
		end
	end
	verbose && println("skipped $border_edges border edges, $flat_edges flat edges")
end

# ╔═╡ 419bc5f1-a68a-4e06-a344-c889bf49ac4f
function vertex_extrusions(ds, ns, mesh, grid, dmax; verbose = true)

	lininds = LinearIndices(size(grid))

	topo = convert(Meshes.HalfEdgeTopology, Meshes.topology(mesh))
	v2e = Meshes.Coboundary{0,1}(topo)
	v2t = Meshes.Coboundary{0,2}(topo)
	e2v = Meshes.Boundary{1,0}(topo)
	t2v = Meshes.Boundary{2,0}(topo)
	covertices(ivs, iv) = ivs[1] == iv ? (ivs[2], ivs[3]) : ivs[2] == iv ? (ivs[1], ivs[3]) : (ivs[1], ivs[2])
	covertex(ivs, iv) = ivs[1] == iv ? ivs[2] : ivs[1]

	vs = Meshes.vertices(mesh)
	plane_points = 0
	convex_points = 0
	concave_points = 0
	saddle_points = 0

	for iv in Meshes.vertices(topo)

		# build pseudo-normal vector
		nα = Meshes.Vec(0,0,0)
		for it in v2t(iv)
			ni = Meshes.normal(mesh[it])
			αi = acos(Meshes.dot((vs[ivt] - vs[iv] for ivt in covertices(t2v(it), iv))...))
			nα += αi * ni
		end
		nα /= Meshes.norm(nα)

		# count points above/below plane of pseudo-normal vector
		nvb, nva = 0, 0
		for ie in v2e(iv)
			loc = Int(sign(Meshes.dot(vs[covertex(e2v(ie), iv)] - vs[iv], nα)))
			nva += max(loc, 0)
			nvb += max(-loc, 0)
		end

		# classify vertices	plane_points = 0
		if nva == nvb == 0
			plane_points += 1
		elseif nva == 0 # convex vertex
			convex_points += 1
		elseif nvb == 0 # concave vertex
			concave_points += 1
		else
			saddle_points += 1
		end

		# define ball for vertex extrusion
		p = vs[iv]
		ball = Meshes.Ball(p, dmax)

		# compute distance for points inside ball
		for i in rasterize(ball, grid)
			x = Meshes.centroid(grid[lininds[i]])
			p2x = x - p
			d = Meshes.norm(p2x) * sign(Meshes.dot(p2x, nα))
			abs(d) <= dmax || continue
			d_old = ds[i]
			if ismissing(d_old) || abs(d) < abs(d_old)
				ds[i] = d
			end
		end
	end
	verbose && println("$plane_points plane points, $saddle_points saddle points, $convex_points convex points, $concave_points concave points")
end

# ╔═╡ e3df4b53-2523-421f-a308-abb4a40565f1
function distancefield(mesh::Meshes.SimpleMesh{Dim,T},
					   grid::Meshes.CartesianGrid{Dim,T},
					   dmax::T; verbose = true) where {Dim,T}
	#ds = Array{Union{Nothing,T}}(nothing, size(grid)...)
	#ns = Array{Union{Nothing,Meshes.Vec{Dim,T}}}(nothing, size(grid)...)
	ds = Array{Union{Missing,T}}(missing, size(grid)...)
	ns = Array{Union{Missing,Meshes.Vec{Dim,T}}}(missing, size(grid)...)
	#ds[3:17, 3:17, 4:17] .= 1

	plane_extrusions(ds, ns, mesh, grid, dmax; verbose = verbose)
	edge_extrusions(ds, ns, mesh, grid, dmax; verbose = verbose)
	vertex_extrusions(ds, ns, mesh, grid, dmax; verbose = verbose)

	Meshes.meshdata(grid, etable=(mesh_distance = ds, mesh_normal = ns))
end

# ╔═╡ 706fec20-d017-4ad8-a51b-9302edf2a01a
function benchmark_sdf(N)
	mesh = loadmesh("single-cube.stl", flip=true)
	grid = Meshes.CartesianGrid((0.,0.,-0.5), (2.,2.,1.5), dims=(N,N,N))
	@benchmark distancefield($mesh, $grid, $0.25, verbose=false)
end

# ╔═╡ 53e05dac-414f-41e5-9b27-c696c17cd9ed
# ╠═╡ disabled = true
#=╠═╡
benchmark_sdf(80)
  ╠═╡ =#

# ╔═╡ c045f033-5ffd-45b1-9156-9b34e1724ac5
let mesh = loadmesh("single-cube.stl", flip=true)
	N = 80
	grid = Meshes.CartesianGrid((0.,0.,-0.5), (2.,2.,1.5), dims=(N,N,N))
	sdf = distancefield(mesh, grid, 0.25)
	d = Meshes.values(sdf).mesh_distance

	iplt = 20
	nx, ny, nz = size(grid)
	x = LinRange(0, 2, 1+2*nx)[2:2:end]
	y = LinRange(0, 2, 1+2*ny)[2:2:end]
	z = LinRange(0, 2, 1+2*nz)[2:2:end]
	fig = Mke.Figure(resolution = (500, 400))
	ax = Mke.Axis(fig[1,1], title="y=$(y[iplt])")
	hm = Mke.heatmap!(x, y, d[:,iplt,:], colorrange = (-0.3,0.3), colormap = :oleron)
	Mke.Colorbar(fig[1, 2], hm)
	fig
end

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
# ╠═9929f2be-277a-11ed-375a-f18a1b52c0ca
# ╟─87aca29f-29df-44ba-9630-cc97e927a339
# ╟─3671d217-c327-431f-9bfd-0f3e7fd4a01a
# ╟─fc1826d8-f6ad-443f-bf65-2d5def8974a7
# ╠═a96a6610-2199-414b-a97b-2f3c4a8d33d2
# ╠═f33325a0-ca9e-4366-bc37-91dc414cad25
# ╠═34b3f98c-1053-4aa7-ae0f-4aff49998c82
# ╟─4ca47d96-e749-4b06-9b1e-8fa2222a0936
# ╠═908217b7-39f8-42bd-870b-92a2e872b7a7
# ╠═a45d21c3-ea3b-4d11-8a4a-be519f1d7719
# ╠═8b9a30c3-f4ab-4367-b2a0-a398da18cd1b
# ╟─6969c27d-6609-4f30-94fd-5d1c9a4fdc91
# ╟─5a1e55bb-7a96-46c4-85d0-cf40606860bf
# ╠═73dfeea2-3faf-458c-ad41-7d592a581479
# ╠═c3ce016a-b6bb-4bf6-aa48-4f5da7818c77
# ╟─82d773d2-a328-41ca-b5c6-1836a393d6da
# ╠═7999f361-6a5b-4712-b7cb-3904fbecd108
# ╟─07d89b37-ffc9-49c1-9f7c-59e26f256fcd
# ╟─17854bbc-543c-42bd-a126-943edda8c75a
# ╟─2d4b6ea3-decb-4f6a-939d-1eafe956828e
# ╟─3772c9e9-2db1-4978-a011-de19799f4e98
# ╟─cb2c81fe-1987-42c9-aec8-badbbd0b1a0a
# ╟─0d7b073f-90a3-4588-b4d1-5cc0b14e2bc5
# ╟─409498df-51ed-49e8-85a9-efcdee401255
# ╠═e0415fec-1d32-458d-ae55-87251486cbcb
# ╠═01843255-f9a3-49a7-a989-b15b2e6932de
# ╟─3c555464-59e2-4c98-bdfe-f64891f6b569
# ╠═bf425675-9ee5-43de-a67a-b9098a751187
# ╠═8148d1c8-8fea-4a7d-a1d4-23a256a10b8a
# ╟─342ed4c3-ffcc-4a75-887d-5c9c6f8318d6
# ╠═1dd29d36-dfd3-4bea-a29c-b190c5713575
# ╟─c49f0f3d-043c-4dc4-aab8-7fa7b5a23dda
# ╠═70f8f6cb-47d2-4254-8032-f93a713f8e74
# ╠═ad748fba-054c-4d74-b589-24dbc0d8afd1
# ╠═8ffe0ca3-07bd-46ff-b525-8890acf98b68
# ╠═a1588120-35c9-4719-8bca-757f56cc288a
# ╠═ba3a5975-631f-44d3-9d7f-96c89927b7a2
# ╟─59335790-2400-4e6b-aff5-5f3837503033
# ╟─86244339-1362-48b0-8342-d4a6528138e3
# ╟─fccdb7f1-4ca3-47ec-85a9-8d6b0b2e23c8
# ╠═e24d37db-79ea-45fc-a5f4-2b372b30f5a8
# ╠═35e13f7a-cde4-430c-9fe8-51c96c2e5af4
# ╠═86e0b953-924e-44d2-acd4-1515e732eb6e
# ╟─9069a1a1-31b2-4e5a-80fd-a410475f4d88
# ╠═3e13d529-2114-4495-a1eb-7669f2a068c7
# ╠═feaefa19-fb04-42c5-bdd7-b095542b9cd2
# ╠═419bc5f1-a68a-4e06-a344-c889bf49ac4f
# ╠═e3df4b53-2523-421f-a308-abb4a40565f1
# ╠═706fec20-d017-4ad8-a51b-9302edf2a01a
# ╠═53e05dac-414f-41e5-9b27-c696c17cd9ed
# ╠═c045f033-5ffd-45b1-9156-9b34e1724ac5
# ╟─ff1d3625-60e0-41dd-b889-0d2606378568
# ╟─a07eeeed-1d9d-4574-bad3-fabc4f0b0af9
# ╟─d4065db9-6b48-4c1e-ae33-ed793fc91c28
# ╠═762e6049-bf06-437c-b990-590a2d025021
# ╠═f809754a-ea9f-4062-8a94-85e64b42c3a2
# ╠═f415eab0-b090-495c-a79c-64a8ad4f3232
# ╟─45c84785-e659-4454-9011-4b3281107498
# ╠═277029cf-135e-465e-9401-f54734caf16e
# ╠═973a989f-8090-4261-af3b-929f4d418633
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
