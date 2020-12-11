# module CyclicRayTracing
module RayTracing

using UnPack
using StaticArrays
using Gridap: Point, norm, num_cells, get_grid
# using GriadpGmsh # ver de usar dev GriadpGmsh para que esto funque?
# using GriadpGmsh.Geometry: UnstructuredDiscreteModel

# IDEAS:
# 1. GridapEmbedded.jl para describir geometrias como open moc, constructive solid geometry.
#    Aca podria terminar implementando como tiene open moc las intersections entre figuras y
#    tracks.
# 2. mirar Meshes.jl y a la organizacion que pertenece. ver funciones de buscar in element, si las tiene
# 3. kdtrees para meshes: https://github.com/KristofferC/NearestNeighbors.jl


include("Quadrature.jl")
include("Track.jl")
include("Mesh.jl")
include("TrackGenerator.jl")

export TrackGenerator
export trace!

end