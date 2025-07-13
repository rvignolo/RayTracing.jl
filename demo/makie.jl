using Plots
using Gridap
using GLMakie
using RayTracing
using DataStructures: CircularBuffer
using Colors
using UnPack

const TAIL_LENGTH = 20000
const BACKGROUND_COLOR = RGBf(0.98, 0.98, 0.98);

# Lightweight mesh drawing function
function draw_lightweight_mesh!(ax, mesh;
    cell_color=RGBf(0.95, 0.95, 0.95),
    edge_color=RGBf(0.8, 0.8, 0.8),
    edge_width=0.5,
    alpha=0.3)
    @unpack model, cell_nodes = mesh
    grid = get_grid(model)
    node_coordinates = Gridap.ReferenceFEs.get_node_coordinates(grid)

    # Pre-allocate arrays to avoid repeated allocations
    cell_count = length(cell_nodes)

    # Pre-allocate arrays for batch operations
    all_polygons = Vector{Vector{Point2f}}(undef, cell_count)

    # Prepare all polygon coordinates in a single pass
    for (cell_id, node_ids) in enumerate(cell_nodes)
        nodes = [node_coordinates[nid] for nid in node_ids]
        all_polygons[cell_id] = Point2f.(getproperty.(nodes, :data))
    end

    # Draw all polygons in a single batch operation
    poly!(ax, all_polygons,
        color=cell_color,
        strokecolor=edge_color,
        strokewidth=edge_width,
        alpha=alpha)
end

# load geometry
jsonfile = joinpath(@__DIR__, "pincell.json")
model = DiscreteModelFromFile(jsonfile)

nφ = 16
δ = 0.08
bc = Reflective
bcs = BoundaryConditions(top=bc, bottom=bc, left=bc, right=bc)
tg = TrackGenerator(model, nφ, δ, bcs=bcs)
trace!(tg)
segmentize!(tg)

fig = Figure(
    resolution=(1000, 1000),
    backgroundcolor=BACKGROUND_COLOR,
    fontsize=16,
    font="Computer Modern",
)

# Create main axis with enhanced styling
ax = Axis(
    fig[1, 1],
    title="Ray Tracing Visualization",
    xlabel="X Position",
    ylabel="Y Position",
    backgroundcolor=:white,
    xgridvisible=false,
    ygridvisible=false,
    xgridcolor=RGBf(0.9, 0.9, 0.9),
    ygridcolor=RGBf(0.9, 0.9, 0.9),
    xgridwidth=1.0,
    ygridwidth=1.0,
    xgridstyle=:dash,
    ygridstyle=:dash,
    xticklabelsize=12,
    yticklabelsize=12,
    xlabelsize=14,
    ylabelsize=14,
    titlesize=18,
    aspect=DataAspect(),
    autolimitaspect=1.0,
    limits=(nothing, nothing, nothing, nothing)
)

# Set axis limits
GLMakie.xlims!(ax, (tg.mesh.bb_min.x, tg.mesh.bb_max.x))
GLMakie.ylims!(ax, (tg.mesh.bb_min.y, tg.mesh.bb_max.y))

# Draw lightweight mesh
draw_lightweight_mesh!(ax, tg.mesh)

# updates the trajectory with a new segment
function update_ray!(trajectory_obs, segment, direction)
    if direction == RayTracing.Forward
        start_point = Point2f(segment.p)
        end_point = Point2f(segment.q)
    else
        start_point = Point2f(segment.q)
        end_point = Point2f(segment.p)
    end

    # add new points to the end of the trajectory
    push!(trajectory_obs[], start_point)
    push!(trajectory_obs[], end_point)
end

# plots a cyclic trajectory with enhanced colors, output can be gif or even mp4
function trajectory(ax, initial_track, initial_direction, output)

    # initialize track and direction
    track = initial_track
    dir = initial_direction

    # a circular buffer maintains its size and pushed values replace the latest one
    x1, y1 = track.p
    trajectory = CircularBuffer{Point2f}(TAIL_LENGTH)
    fill!(trajectory, Point2f(x1, y1))
    trajectory_obs = Observable(trajectory)

    # Draw the trajectory with enhanced styling
    lines!(ax, trajectory_obs;
        linewidth=2.5,
        color=to_color(:black),
        linestyle=:solid)

    record(fig, output) do io

        # Only update observable every N segments, for performance reasons
        update_counter = 0
        update_frequency = 5

        while true

            # the segments are stored in reverse order for backward tracks
            segments = dir == RayTracing.Backward ? reverse(track.segments) : track.segments

            # for each segment, update the trajectory and record the frame
            for (i, segment) in enumerate(segments)

                update_ray!(trajectory_obs, segment, dir)
                update_counter += 1

                # Only update the observable periodically to reduce overhead
                if update_counter % update_frequency == 0
                    trajectory_obs[] = trajectory_obs[]
                end

                # record all segments, but we could record only if a condition was met (e.g every 10 segments)
                recordframe!(io)
            end

            # Force final update
            trajectory_obs[] = trajectory_obs[]

            # to separate tracks into segments
            push!(trajectory_obs[], Point2f(NaN, NaN))

            # update the track and direction
            if dir == RayTracing.Forward
                dir = RayTracing.dir_next_track_fwd(track)
                track = track.next_track_fwd
            else
                dir = RayTracing.dir_next_track_bwd(track)
                track = track.next_track_bwd
            end

            # stop the animation when the track returns to the initial track
            track.uid != initial_track.uid || break
        end
    end
end

# Plot a single track, purposefully selected to look nice
initial_track = tg.tracks[2][1]
trajectory(ax, initial_track, RayTracing.Forward, "cyclic_track_with_mesh.gif")

# This also works and would plot all tracks simultaneously, but it has limited performance
# for (i, azimuthal_tracks) in enumerate(tg.tracks)
#     @async trajectory(ax, first(azimuthal_tracks), RayTracing.Forward, "cyclic_track.gif")
# end
