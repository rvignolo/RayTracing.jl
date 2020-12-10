# module CyclicRayTracing
module RayTracing

using UnPack
using StaticArrays
using Gridap: Point


# IDEAS:
# 1. GridapEmbedded.jl para describir geometrias como open moc, constructive solid geometry.
#    Aca podria terminar implementando como tiene open moc las intersections entre figuras y
#    tracks.
# 2. mirar Meshes.jl y a la organizacion que pertenece.
# 3. anotar todas las funciones que necesito, como mesh_find_element(mesh, x) de wasora
#    (encontrar la cell que tiene el point), etc, y consultar si existen.


include("Quadrature.jl")
include("Track.jl")
include("TrackGenerator.jl")

end