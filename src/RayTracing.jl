module RayTracing

using UnPack
using RecipesBase
using StaticArrays
using IntervalSets
using LinearAlgebra
using NearestNeighbors
using Gridap: VectorValue, num_cells, get_grid
using Gridap.Geometry: UnstructuredDiscreteModel, UnstructuredGrid, get_grid_topology

using Gridap.ReferenceFEs: get_faces, get_node_coordinates, num_cell_dims

import Base: show
import NearestNeighbors: KDTree

include("quadrature.jl")
include("track.jl")
include("mesh.jl")
include("tracer.jl")
include("plot_recipes.jl")

export TrackGenerator
export trace!
export segmentize!

end
