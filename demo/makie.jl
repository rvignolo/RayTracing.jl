using Gridap
using GLMakie
using RayTracing
using DataStructures: CircularBuffer
using Colors
using UnPack

const TAIL_LENGTH = 20000
const BACKGROUND_COLOR = RGBf(0.98, 0.98, 0.98);

# Transport problem.
jsonfile = joinpath(@__DIR__, "pincell.json")
model = DiscreteModelFromFile(jsonfile)
nφ = 16
δ = 0.08
bc = Reflective
bcs = BoundaryConditions(top=bc, bottom=bc, left=bc, right=bc)
tg = TrackGenerator(model, nφ, δ, bcs=bcs)
trace!(tg)
segmentize!(tg)

# Draw a lightweight mesh.
function draw_lightweight_mesh!(ax, mesh;
    cell_color=RGBf(0.95, 0.95, 0.95),
    edge_color=RGBf(0.8, 0.8, 0.8),
    edge_width=0.5,
    alpha=0.3)
    @unpack model, cell_nodes = mesh
    grid = get_grid(model)
    node_coordinates = Gridap.ReferenceFEs.get_node_coordinates(grid)

    face_labeling = get_face_labeling(model)
    cell_tags = Gridap.Geometry.get_face_tag(face_labeling, 2)

    # Preallocate arrays to avoid repeated allocations.
    cell_count = length(cell_nodes)
    cell_colors = Vector{RGBf}(undef, cell_count)
    all_polygons = Vector{Vector{Point2f}}(undef, cell_count)

    # Prepare all polygon coordinates in a single pass.
    for (cell_id, node_ids) in enumerate(cell_nodes)
        nodes = [node_coordinates[nid] for nid in node_ids]
        all_polygons[cell_id] = Point2f.(getproperty.(nodes, :data))

        # Color by cell type, darkening higher tags.
        cell_tag = cell_tags[cell_id]
        cell_color′ = cell_color * (cell_tag / maximum(cell_tags))
        cell_colors[cell_id] = cell_color′
    end

    # Draw all polygons in a single batch operation.
    poly!(ax, all_polygons,
        color=cell_colors,
        strokecolor=edge_color,
        strokewidth=edge_width,
        alpha=alpha)
end

# Update the trajectory with a new segment.
function update_ray!(trajectory_obs, segment, direction)
    if direction == RayTracing.Forward
        start_point = Point2f(segment.p)
        end_point = Point2f(segment.q)
    else
        start_point = Point2f(segment.q)
        end_point = Point2f(segment.p)
    end

    # Add new points to the end of the trajectory.
    push!(trajectory_obs[], start_point)
    push!(trajectory_obs[], end_point)
end

# Plot a cyclic trajectory; output can be a GIF or an MP4.
function trajectory(fig, ax, initial_track, initial_direction, output)

    # Initialize the track and direction.
    track = initial_track
    dir = initial_direction

    # A circular buffer keeps a fixed tail length as new values are pushed.
    x1, y1 = track.p
    trajectory = CircularBuffer{Point2f}(TAIL_LENGTH)
    fill!(trajectory, Point2f(x1, y1))
    trajectory_obs = Observable(trajectory)

    # Draw the trajectory.
    lines!(ax, trajectory_obs;
        linewidth=2.5,
        color=to_color(:black),
        linestyle=:solid)

    record(fig, output, framerate=60) do io

        # Only update the observable every N segments for performance.
        update_counter = 0
        update_frequency = 5

        while true

            # Segments are stored in reverse order for backward tracks.
            segments = dir == RayTracing.Backward ? reverse(track.segments) : track.segments

            # Update the trajectory and record one frame per segment.
            for (i, segment) in enumerate(segments)

                update_ray!(trajectory_obs, segment, dir)
                update_counter += 1

                # Update the observable periodically to reduce overhead.
                if update_counter % update_frequency == 0
                    trajectory_obs[] = trajectory_obs[]
                end

                # Record every segment; this could be throttled, for example every 10 segments.
                recordframe!(io)
            end

            # Force a final update.
            trajectory_obs[] = trajectory_obs[]

            # Separate tracks in the rendered trajectory.
            push!(trajectory_obs[], Point2f(NaN, NaN))

            # Update the track and direction.
            if dir == RayTracing.Forward
                dir = RayTracing.dir_next_track_fwd(track)
                track = track.next_track_fwd
            else
                dir = RayTracing.dir_next_track_bwd(track)
                track = track.next_track_bwd
            end

            # Stop when the trajectory returns to the initial track.
            track.uid != initial_track.uid || break
        end
    end
end

fig = Figure(
    resolution=(1000, 1000),
    backgroundcolor=BACKGROUND_COLOR,
    fontsize=16,
    font="Computer Modern",
)

# Create the main axis.
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

# Set axis limits.
GLMakie.xlims!(ax, (tg.mesh.bb_min.x, tg.mesh.bb_max.x))
GLMakie.ylims!(ax, (tg.mesh.bb_min.y, tg.mesh.bb_max.y))

# Draw the lightweight mesh.
draw_lightweight_mesh!(ax, tg.mesh)

# Plot one track selected for a clear visualization.
initial_track = tg.tracks[2][1]
trajectory(fig, ax, initial_track, RayTracing.Forward, joinpath(@__DIR__, "cyclic_track_with_mesh.gif"))

# This also plots all tracks simultaneously, but performance is poor.
# for (i, azimuthal_tracks) in enumerate(tg.tracks)
#     @async trajectory(ax, first(azimuthal_tracks), RayTracing.Forward, "cyclic_track.gif")
# end
