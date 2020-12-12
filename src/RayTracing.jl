module RayTracing

using UnPack
using RecipesBase
using StaticArrays
using Gridap: Point, norm, num_cells, get_grid
# using GriadpGmsh # ver de usar dev GriadpGmsh para que esto funque?
# using GriadpGmsh.Geometry: UnstructuredDiscreteModel

include("Quadrature.jl")
include("Track.jl")
include("Mesh.jl")
include("TrackGenerator.jl")

export TrackGenerator
export trace!

end
