module RayTracing

using UnPack
using RecipesBase
using StaticArrays
using IntervalSets
using NearestNeighbors
using Gridap: Point, norm, num_cells, get_grid
using Gridap.Geometry: UnstructuredGrid, get_grid_topology
# using GriadpGmsh # ver de usar dev GriadpGmsh para que esto funque?

using Gridap.ReferenceFEs: get_faces, get_node_coordinates, num_cell_dims

import Base: show
import NearestNeighbors: KDTree

include("quadrature.jl")
include("track.jl")
include("mesh.jl")
include("tracer.jl")

export TrackGenerator
export trace!

end
