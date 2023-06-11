using Plots
using Gridap
using GLMakie
using RayTracing
using DataStructures: CircularBuffer


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

# proceed to segmentation
segmentize!(tg)

# plot
nφ = 16
δ = 0.08
bc = Reflective
bcs = BoundaryConditions(top=bc, bottom=bc, left=bc, right=bc)
tg = TrackGenerator(model, nφ, δ, bcs=bcs)
trace!(tg)
segmentize!(tg)

fig = Figure();
display(fig);
ax = Axis(fig[1, 1])
GLMakie.xlims!(ax, (tg.mesh.bbmin.x, tg.mesh.bbmax.x))
GLMakie.ylims!(ax, (tg.mesh.bbmin.y, tg.mesh.bbmax.y))

# this plots a cyclic trajectory. I have selected a track that looks nice
trajectory(ax, tg.tracks[2][1], RayTracing.Forward, "cyclic_track.gif")

# this also works and would plot all tracks simultaneously
# for azimutal_tracks in tg.tracks
#     @async trajectory(ax, azimutal_tracks[1], RayTracing.Forward, "cyclic_track.mp4")
# end

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

# plots a ciclyc trajectory
function trajectory(ax, initial_track, initial_dir, output)

    track = initial_track
    dir = initial_dir

    x1, y1 = track.p
    tail = 100_000
    traj = CircularBuffer{Point2f}(tail)
    fill!(traj, Point2f(x1, y1))
    traj = Observable(traj)

    c = to_color(:black)
    tailcol = [RGBAf(c.r, c.g, c.b, (i / tail)^2) for i in 1:tail]
    lines!(ax, traj; linewidth=1, color=tailcol)
    # lines!(ax, traj; linestyle = :dot, linewidth = 3, color = tailcol)

    # changing the extension here to mp4 would produce a video
    record(fig, output) do io

        while true
            segments = dir == RayTracing.Backward ? reverse(track.segments) : track.segments
            for segment in segments
                animstep!(segment, dir, traj)
                # sleep(0.002)
                recordframe!(io) # record all segments, but we could record only if a condition was met
            end

            # to separate tracks into segments
            push!(traj[], Point2f(NaN, NaN))

            if dir == RayTracing.Forward
                dir = RayTracing.dir_next_track_fwd(track)
                track = track.next_track_fwd
            else
                dir = RayTracing.dir_next_track_bwd(track)
                track = track.next_track_bwd
            end

            track.uid != initial_track.uid || break
        end
    end
end
