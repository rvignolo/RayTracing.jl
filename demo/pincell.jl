using Plots
using RayTracing
using GridapGmsh: GmshDiscreteModel

mshfile = joinpath(@__DIR__,"pincell.msh")
model = GmshDiscreteModel(mshfile; renumber=true)

# number of azimuthal angles
nφ = 8

# azimuthal spacing
δ = 0.04

# initialize track generator
tg = TrackGenerator(model, nφ, δ)

# perform ray tracing
trace!(tg)

# plot tracks
plot(tg, dpi=300)

savefig("pincell.pdf")
savefig("pincell.png")