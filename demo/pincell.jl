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
function trajectory(ax, initial_track, initial_dir)

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

    record(fig, "cyclic_track.gif") do io

        while true
            segments = dir == RayTracing.Backward ? reverse(track.segments) : track.segments
            j = 0
            for segment in segments
                animstep!(segment, dir, traj)
                sleep(0.005)
                j += 1
                # if j % 2 == 0
                recordframe!(io)
                # end
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

            track.uid != initial_track.uid || break
        end

    end

end

trajectory(ax, tg.tracks[2][1], RayTracing.Forward)

# elegir varios tracks que tengan diferentes azimuthal
for azimutal_tracks in tg.tracks[2]
    @async trajectory(ax, azimutal_tracks[1], RayTracing.Forward)
end


# record(fig, "video.mp4", frames; framerate = 60) do i # i = frame number
#     for j in 1:5 # step 5 times per frame
#         animstep!(segment, dir,traj)
#     end
#     # any other manipulation of the figure here...
# end



# mesh plot
# @recipe function plot(mesh::Mesh)
mesh1 = tg.mesh
using UnPack
@unpack cell_nodes, model = mesh1

grid = get_grid(model)
nodes = RayTracing.get_node_coordinates(grid)

# number of cells and number of nodes per cell
nc = length(cell_nodes)
nn = all(l -> isequal(l, length(cell_nodes[1])), length.(cell_nodes)) ? length(cell_nodes[1]) : error("error")

x = Vector{Float64}(undef, (nn + 1) * nc + nc + nc) # TODO: compute Float64
y = Vector{Float64}(undef, (nn + 1) * nc + nc + nc)

k = 1
for (i, nodes_ids) in enumerate(cell_nodes)

    for (j, node_id) in enumerate(nodes_ids)

        x[k] = nodes[node_id][1]
        y[k] = nodes[node_id][2]

        k += 1
    end

    x[k] = nodes[nodes_ids[1]][1]
    y[k] = nodes[nodes_ids[1]][2]

    k += 1

    x[k] = NaN
    y[k] = NaN

    k += 1
end

lines!(ax, x, y, linewidth=0.3, color=:black)

#     linecolor   --> :black
#     seriestype  :=  :path
#     linewidth   --> 0.2
#     legend      --> false
#     framestyle  :=  :none

#     return (x, y)
# end