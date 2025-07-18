module RayTracing

using UnPack
using RecipesBase
using StaticArrays
using IntervalSets
using LinearAlgebra
using NearestNeighbors
using Gridap: VectorValue, num_cells, get_grid
using Gridap.ReferenceFEs: get_faces, get_node_coordinates, num_cell_dims
using Gridap.Geometry:
    UnstructuredDiscreteModel,
    UnstructuredGrid,
    get_grid_topology,
    get_cell_node_ids

import Base: show, in
import NearestNeighbors: KDTree
import Gridap.Geometry: num_cells
import Gridap.ReferenceFEs: num_dims, num_nodes

include("point.jl")
include("segment.jl")
include("boundary.jl")
include("track.jl")
include("mesh.jl")
include("intersection.jl")
include("azimuthal_quad.jl")
include("trackgenerator.jl")
include("plot_recipes.jl")

export BoundaryConditions, Vacuum, Reflective, Periodic
export TrackGenerator
export trace!
export segmentize!

end
