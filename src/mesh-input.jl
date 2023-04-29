module MeshInput

using MeshIO: MeshIO
using FileIO: FileIO
using Meshes: Point, SimpleMesh, connect

export loadmesh

# helper to load various mesh files to Meshes.jl format
function loadmesh(path; flip = false)
    mesh = FileIO.load(path)
    vertices = collect(Set(Point{3,Float64}(v) for v in MeshIO.GeometryBasics.coordinates(mesh)))
    indices = Dict(p => i for (i, p) in enumerate(vertices))
    connectivities = map(mesh) do el
        i, j, k = (indices[Point{3,Float64}(p)] for p in el)
        connect(flip ? (i,k,j) : (i,j,k))
    end

    SimpleMesh(vertices, connectivities, relations=true)
end

end # module MeshInput
