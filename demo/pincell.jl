using Plots
using RayTracing
using Gridap

# load geometry
jsonfile = joinpath(@__DIR__,"pincell.json")
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
# plot(tg, dpi=300, size=(250,250), palette=:Paired_4, background_color=:transparent)

# savefig("pincell.svg")
# savefig("pincell.pdf")
# savefig("pincell.png")

# proceed to segmentation
segmentize!(tg)

# plot mesh
# plot(tg.mesh)

# plot segments
# plot(tg.tracks_by_uid)

# savefig("segments.svg")










using GLMakie
using DataStructures: CircularBuffer

nφ = 8
δ = 0.1
bc = Reflective
bcs = BoundaryConditions(top=bc, bottom=bc, left=bc, right=bc)
tg = TrackGenerator(model, nφ, δ, bcs=bcs)
trace!(tg)
segmentize!(tg)

fig = Figure(); display(fig)
ax = Axis(fig[1,1])
GLMakie.xlims!(ax, (tg.mesh.bbmin.x, tg.mesh.bbmax.x))
GLMakie.ylims!(ax, (tg.mesh.bbmin.y, tg.mesh.bbmax.y))

track = first(tg.tracks_by_uid)
x1, y1 = track.p
tail = 50_000
traj = CircularBuffer{Point2f}(tail)
fill!(traj, Point2f(x1, y1))
traj = Observable(traj)

c = to_color(:black);
tailcol = [RGBAf(c.r, c.g, c.b, (i/tail)^2) for i in 1:tail];
lines!(ax, traj; linewidth = 1, color = tailcol)
# lines!(ax, traj; linestyle = :dot, linewidth = 3, color = tailcol)

function animstep!(segment, dir, traj)
    if dir == RayTracing.Forward
        px, py = segment.p
        qx, qy = segment.q
    else
        px, py = segment.q
        qx, qy = segment.p
    end
    push!(traj[], Point2f(px, py))
    push!(traj[], Point2f(qx, qy))
    traj[] = traj[]
end


# esta es una funcion que plotea una trajectoria que es cyclica. creo que deberia lanzar una
# trajectoria por n_track_x y n_track_y
track = first(tg.tracks_by_uid)
dir = RayTracing.Forward
while(track !== nothing)

    segments = dir == RayTracing.Backward ? reverse(track.segments) : track.segments
    for segment in segments
        animstep!(segment, dir, traj)
        # sleep(0.02)
        sleep(0.001)
    end

    # to separate tracks in segments
    push!(traj[], Point2f(NaN, NaN))

    # @show dir
    # @show track.track_idx

    if dir == RayTracing.Forward
        dir = RayTracing.dir_next_track_fwd(track)
        track = track.next_track_fwd
    else
        dir = RayTracing.dir_next_track_bwd(track)
        track = track.next_track_bwd
    end

end