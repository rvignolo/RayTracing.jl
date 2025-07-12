using Plots
using Gridap
using GLMakie
using RayTracing
using DataStructures: CircularBuffer
using Colors
using LinearAlgebra
using UnPack

# Enhanced Configuration
const FPS = 30
const TAIL_LENGTH = 300
const RAY_WIDTH = 3.0
const MESH_ALPHA = 0.2
const BACKGROUND_COLOR = RGBf(0.98, 0.98, 0.98)
const MESH_COLOR = RGBf(0.8, 0.8, 0.8)
const BOUNDARY_COLOR = RGBf(0.3, 0.3, 0.3)
const RAY_COLORS = [
    RGBf(0.9, 0.2, 0.2),  # Bright Red
    RGBf(0.2, 0.7, 0.9),  # Bright Blue
    RGBf(0.2, 0.9, 0.2),  # Bright Green
    RGBf(0.9, 0.7, 0.2),  # Bright Orange
    RGBf(0.7, 0.2, 0.9),  # Bright Purple
    RGBf(0.9, 0.2, 0.7),  # Bright Pink
]

# load geometry
jsonfile = joinpath(@__DIR__, "pincell.json")
model = DiscreteModelFromFile(jsonfile)

# Setup ray tracing
nφ = 8
δ = 1
bc = Reflective
bcs = BoundaryConditions(top=bc, bottom=bc, left=bc, right=bc)

tg = TrackGenerator(model, nφ, δ, bcs=bcs)
trace!(tg)
segmentize!(tg)

# Create enhanced figure with multiple panels
fig = Figure(
    resolution=(1600, 1000),
    backgroundcolor=BACKGROUND_COLOR,
    fontsize=16
)

# Main visualization panel
ax_main = Axis(
    fig[1, 1],
    title="Ray Tracing Visualization with Mesh",
    xlabel="X Position",
    ylabel="Y Position",
    backgroundcolor=:white,
    xgridvisible=false,
    ygridvisible=false,
    aspect=DataAspect()
)

# Statistics panel
ax_stats = Axis(
    fig[1, 2],
    title="Ray Statistics",
    xlabel="Time Steps",
    ylabel="Distance Traveled",
    backgroundcolor=:white
)

# Set main axis limits
Makie.xlims!(ax_main, (tg.mesh.bb_min.x, tg.mesh.bb_max.x))
Makie.ylims!(ax_main, (tg.mesh.bb_min.y, tg.mesh.bb_max.y))

# Enhanced mesh drawing with boundary highlighting
function draw_enhanced_mesh!(ax, mesh)
    @unpack model, cell_nodes = mesh
    grid = get_grid(model)
    node_coordinates = Gridap.ReferenceFEs.get_node_coordinates(grid)

    # Draw cells
    for (cell_id, node_ids) in enumerate(cell_nodes)
        nodes = [node_coordinates[nid] for nid in node_ids]
        poly_x = [node[1] for node in nodes]
        poly_y = [node[2] for node in nodes]

        poly!(ax, Point2f.(poly_x, poly_y),
            color=MESH_COLOR,
            strokecolor=BOUNDARY_COLOR,
            strokewidth=1.0,
            alpha=MESH_ALPHA)
    end

    # Highlight domain boundaries
    bb_min, bb_max = mesh.bb_min, mesh.bb_max
    boundary_points = [
        Point2f(bb_min.x, bb_min.y), Point2f(bb_max.x, bb_min.y),
        Point2f(bb_max.x, bb_max.y), Point2f(bb_min.x, bb_max.y),
        Point2f(bb_min.x, bb_min.y)
    ]

    lines!(ax, boundary_points,
        color=BOUNDARY_COLOR,
        linewidth=3.0,
        linestyle=:solid)
end

# Draw enhanced mesh
draw_enhanced_mesh!(ax_main, tg.mesh)

# Function to create enhanced animated ray
function create_enhanced_ray(ax, track, direction, color, ray_id)
    # Create trajectory buffer
    traj = CircularBuffer{Point2f}(TAIL_LENGTH)
    fill!(traj, Point2f(track.p))
    traj_obs = Observable(traj)

    # Create enhanced color gradient
    trail_colors = [RGBAf(color.r, color.g, color.b, (i / TAIL_LENGTH)^2) for i in 1:TAIL_LENGTH]

    # Draw the trail with glow effect
    lines!(ax, traj_obs,
        linewidth=RAY_WIDTH + 2,
        color=RGBAf(color.r, color.g, color.b, 0.3),
        linestyle=:solid)
    lines!(ax, traj_obs,
        linewidth=RAY_WIDTH,
        color=trail_colors,
        linestyle=:solid)

    # Add animated ray head with glow
    ray_head = Observable(Point2f(track.p))
    Makie.scatter!(ax, ray_head,
        color=RGBAf(color.r, color.g, color.b, 0.6),
        markersize=12,
        marker=:circle)
    Makie.scatter!(ax, ray_head,
        color=color,
        markersize=8,
        marker=:circle,
        strokecolor=:white,
        strokewidth=2)

    return traj_obs, ray_head
end

# Function to update ray with statistics
function update_enhanced_ray!(traj_obs, ray_head, segment, direction, stats_obs)
    if direction == RayTracing.Forward
        start_point = Point2f(segment.p)
        end_point = Point2f(segment.q)
    else
        start_point = Point2f(segment.q)
        end_point = Point2f(segment.p)
    end

    push!(traj_obs[], start_point)
    push!(traj_obs[], end_point)
    traj_obs[] = traj_obs[]
    ray_head[] = end_point

    # Update statistics
    segment_length = norm(end_point - start_point)
    push!(stats_obs[], segment_length)
    stats_obs[] = stats_obs[]
end

# Function to create comprehensive animation
function create_comprehensive_animation(ax_main, ax_stats, tracks, directions, output_file)
    # Create multiple rays
    rays = []
    stats_observables = []

    for (i, (track, dir)) in enumerate(zip(tracks, directions))
        color = RAY_COLORS[mod(i - 1, length(RAY_COLORS))+1]
        traj_obs, ray_head = create_enhanced_ray(ax_main, track, dir, color, i)

        # Create statistics observable
        stats = Observable(Float32[])
        push!(stats_observables, stats)

        # Plot statistics line
        lines!(ax_stats, stats,
            color=color,
            linewidth=2,
            label="Ray $i")

        push!(rays, (track, dir, traj_obs, ray_head, color, stats))
    end

    # Add legend to stats
    axislegend(ax_stats, position=:rt)

    # Animation loop
    record(fig, output_file, framerate=FPS) do io
        max_iterations = 800
        iteration = 0

        while iteration < max_iterations
            iteration += 1

            # Update each ray
            for (track, dir, traj_obs, ray_head, color, stats) in rays
                segments = dir == RayTracing.Backward ? reverse(track.segments) : track.segments

                for segment in segments
                    update_enhanced_ray!(traj_obs, ray_head, segment, dir, stats)
                    recordframe!(io)
                end

                # Add gap between tracks
                push!(traj_obs[], Point2f(NaN, NaN))

                # Move to next track
                if dir == RayTracing.Forward
                    dir = RayTracing.dir_next_track_fwd(track)
                    track = track.next_track_fwd
                else
                    dir = RayTracing.dir_next_track_bwd(track)
                    track = track.next_track_bwd
                end

                # Update ray data
                track_idx = findfirst(r -> r[1] == track, rays)
                if track_idx !== nothing
                    rays[track_idx] = (track, dir, traj_obs, ray_head, color, stats)
                end
            end

            # Check if all rays have completed a cycle
            if all(r -> r[1].uid == rays[1][1].uid, rays)
                break
            end
        end
    end
end

# Add information text
info_text = """
Ray Tracing Parameters:
• Azimuthal angles: $nφ
• Spacing: $δ
• Total tracks: $(tg.n_total_tracks)
• Boundary: Reflective
"""

text!(fig[2, 1], info_text,
    position=(0, 0),
    textsize=12,
    color=:black,
    align=(:left, :top))

# Create comprehensive animation
println("Creating comprehensive ray tracing animation...")
tracks_to_animate = [tg.tracks[i][1] for i in 1:min(3, length(tg.tracks))]
directions = [RayTracing.Forward for _ in tracks_to_animate]

create_comprehensive_animation(ax_main, ax_stats, tracks_to_animate, directions, "comprehensive_ray_tracing.gif")

println("Enhanced animation complete! Check 'comprehensive_ray_tracing.gif'")