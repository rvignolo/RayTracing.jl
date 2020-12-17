using Plots
using RayTracing
using GridapGmsh: GmshDiscreteModel

mshfile = joinpath(@__DIR__,"pincell.msh")
model = GmshDiscreteModel(mshfile; renumber=true)

# number of azimuthal angles
nφ = 8

# azimuthal spacing
δ = 0.02

# initialize track generator
tg = TrackGenerator(model, nφ, δ)

# perform ray tracing
trace!(tg)

# plot tracks
# plot(tg, dpi=300, size=(250,250), linecolor=:turquoise, background_color=:transparent)
plot(tg, dpi=300, size=(250,250), palette=:Paired_4, background_color=:transparent)

# savefig("pincell.svg")
# savefig("pincell.pdf")
savefig("pincell.png")

# proceed to segmentation
segmentize!(tg)

# plot mesh
plot(tg.mesh)

# plot segments
for i in 1:tg.n_total_tracks; plot!(tg.tracks_by_uid[i].segments); end;

savefig("mesh-segments.svg")