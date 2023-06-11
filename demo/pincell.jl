using Plots
using RayTracing
using Gridap

# load geometry
jsonfile = joinpath(@__DIR__, "pincell.json")
model = DiscreteModelFromFile(jsonfile)

# number of azimuthal angles
nφ = 8

# azimuthal spacing
δ = 0.02

# boundary conditions
bcs = BoundaryConditions(top=Reflective, bottom=Reflective, left=Reflective, right=Reflective)

# initialize track generator
tg = TrackGenerator(model, nφ, δ, bcs=bcs)

# perform ray tracing
trace!(tg)

# plot tracks
# plot(tg, dpi=300, size=(250,250), linecolor=:turquoise, background_color=:transparent)
plot(tg, dpi=300, size=(250, 250), palette=:Paired_4, background_color=:transparent)

# savefig("pincell.svg")
# savefig("pincell.pdf")
# savefig("pincell.png")

# proceed to segmentation
segmentize!(tg)

# plot mesh
plot(tg.mesh)

# plot segments
plot(tg.tracks_by_uid)

# savefig("segments.svg")