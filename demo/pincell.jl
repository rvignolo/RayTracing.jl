using Gridap
using Plots
using RayTracing

const OUTDIR = @__DIR__
const IMAGE_SIZE = (900, 900)
const GIF_SIZE = (640, 640)
const GIF_FRAMES = 96
const GIF_FPS = 8
const MATERIAL_COLORS = ["#E76F51", "#F2C14E", "#2A9D8F"]

gr()
default(
    dpi=240,
    fontfamily="Computer Modern",
    background_color=:white,
    foreground_color_subplot="#2F3437",
    grid=false,
    legend=false,
)

model = DiscreteModelFromFile(joinpath(OUTDIR, "pincell.json"))
tg = TrackGenerator(model, 8, 0.08; bcs=reflective_boundaries())
trace!(tg)
segmentize!(tg)

function cell_tags(mesh)
    labeling = get_face_labeling(mesh.model)
    return Gridap.Geometry.get_face_tag(labeling, 2)
end

function mesh_limits(mesh)
    return (
        mesh.bb_min[1], mesh.bb_max[1],
        mesh.bb_min[2], mesh.bb_max[2],
    )
end

function polish!(plt, mesh; size=IMAGE_SIZE)
    xmin, xmax, ymin, ymax = mesh_limits(mesh)
    margin = 0.025 * max(xmax - xmin, ymax - ymin)
    plot!(
        plt;
        aspect_ratio=:equal,
        framestyle=:none,
        ticks=false,
        xlims=(xmin - margin, xmax + margin),
        ylims=(ymin - margin, ymax + margin),
        size=size,
        title="",
        xlabel="",
        ylabel="",
    )
    return plt
end

function plot_material_mesh(mesh; fillalpha=0.9, edgealpha=0.55, size=IMAGE_SIZE)
    nodes = Gridap.ReferenceFEs.get_node_coordinates(get_grid(mesh.model))
    tags = cell_tags(mesh)
    plt = plot()

    for (node_ids, tag) in zip(mesh.cell_nodes, tags)
        ids = RayTracing.ordered_node_ids(mesh, node_ids)
        xs = [nodes[id][1] for id in ids]
        ys = [nodes[id][2] for id in ids]
        color = MATERIAL_COLORS[mod1(tag, length(MATERIAL_COLORS))]
        plot!(
            plt,
            Shape(xs, ys);
            fillcolor=color,
            fillalpha=fillalpha,
            linecolor="#FFFFFF",
            linealpha=edgealpha,
            linewidth=0.25,
        )
    end

    return polish!(plt, mesh; size=size)
end

function plot_material_regions(mesh; size=GIF_SIZE)
    plt = plot_material_mesh(mesh; fillalpha=0.30, edgealpha=0.0, size=size)
    return polish!(plt, mesh; size=size)
end

function plot_material_context(mesh; size=GIF_SIZE)
    plt = plot_material_mesh(mesh; fillalpha=0.30, edgealpha=0.0, size=size)
    plot!(
        plt,
        mesh;
        linecolor="#98A2AE",
        linealpha=0.30,
        linewidth=0.16,
    )
    return polish!(plt, mesh; size=size)
end

function save_static_assets(tg)
    geometry_plot = plot_material_regions(tg.mesh)
    savefig(geometry_plot, joinpath(OUTDIR, "pincell-geometry.png"))

    mesh_plot = plot_material_mesh(tg.mesh)
    savefig(mesh_plot, joinpath(OUTDIR, "pincell-msh.png"))

    tracks_plot = plot_material_mesh(tg.mesh; fillalpha=0.18, edgealpha=0.25)
    plot!(
        tracks_plot,
        tg;
        linecolor="#1F2933",
        linealpha=0.46,
        linewidth=0.42,
    )
    polish!(tracks_plot, tg.mesh)
    savefig(tracks_plot, joinpath(OUTDIR, "pincell-tracks.png"))

    segments_plot = plot(tg.mesh; linecolor="#D6DBDF", linealpha=0.65, linewidth=0.18)
    plot!(
        segments_plot,
        tg.tracks_by_uid;
        color=:turbo,
        colorbar=false,
        linewidth=0.58,
        linealpha=0.9,
    )
    polish!(segments_plot, tg.mesh)
    savefig(segments_plot, joinpath(OUTDIR, "pincell-segments.png"))
end

function cyclic_path(initial_track, initial_direction)
    xs = Float64[]
    ys = Float64[]

    track = initial_track
    direction = initial_direction

    for _ in 1:10_000
        segments = direction == RayTracing.Backward ? reverse(track.segments) : track.segments
        for segment in segments
            p, q = direction == RayTracing.Forward ? (segment.p, segment.q) : (segment.q, segment.p)
            if isempty(xs)
                push!(xs, p[1])
                push!(ys, p[2])
            end
            push!(xs, q[1])
            push!(ys, q[2])
        end

        if direction == RayTracing.Forward
            direction = RayTracing.dir_next_track_fwd(track)
            track = track.next_track_fwd
        else
            direction = RayTracing.dir_next_track_bwd(track)
            track = track.next_track_bwd
        end

        track.uid == initial_track.uid && break
    end

    return xs, ys
end

function plot_ray_tail!(plt, xs, ys, first_idx, last_idx)
    first_idx == last_idx && return plt

    plot!(
        plt,
        xs[first_idx:last_idx],
        ys[first_idx:last_idx];
        seriescolor=:white,
        linecolor=:white,
        linewidth=2.6,
    )

    return plt
end

function save_cycle_animation(tg, output)
    xs, ys = cyclic_path(tg.tracks[2][1], RayTracing.Forward)
    nframes = GIF_FRAMES
    tail = min(260, length(xs))
    base = plot_material_context(tg.mesh)

    animation = @animate for frame in 1:nframes
        idx = max(2, round(Int, 1 + (length(xs) - 1) * (frame - 1) / (nframes - 1)))
        tail_idx = max(1, idx - tail)

        plt = deepcopy(base)
        plot_ray_tail!(plt, xs, ys, tail_idx, idx)
        scatter!(
            plt,
            [xs[idx]],
            [ys[idx]];
            markercolor=:white,
            markerstrokecolor="#334155",
            markerstrokewidth=0.8,
            markersize=4.0,
        )
    end

    gif(animation, joinpath(OUTDIR, output); fps=GIF_FPS)
end

save_static_assets(tg)
save_cycle_animation(tg, "cyclic_track_with_mesh.gif")
