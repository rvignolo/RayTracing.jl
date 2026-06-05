pushfirst!(LOAD_PATH, normpath(joinpath(@__DIR__, "..")))

using Gridap
using Plots
using RayTracing

const OUTDIR = @__DIR__
const IMAGE_SIZE = (900, 900)
const GIF_SIZE = (640, 640)
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

function plot_mesh_outline(mesh; size=GIF_SIZE)
    plt = plot(
        mesh;
        linecolor="#D1D7DC",
        linealpha=0.58,
        linewidth=0.18,
    )
    return polish!(plt, mesh; size=size)
end

function plot_domain(mesh; size=GIF_SIZE)
    xmin, xmax, ymin, ymax = mesh_limits(mesh)
    plt = plot()
    plot!(
        plt,
        Shape([xmin, xmax, xmax, xmin], [ymin, ymin, ymax, ymax]);
        fillcolor=:white,
        fillalpha=1.0,
        linecolor="#C5CDD3",
        linewidth=1.4,
    )
    return polish!(plt, mesh; size=size)
end

function save_static_assets(tg)
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

function save_cycle_animation(tg, output; with_mesh)
    xs, ys = cyclic_path(tg.tracks[2][1], RayTracing.Forward)
    nframes = 72
    tail = min(280, length(xs))
    base_plot = with_mesh ? plot_mesh_outline : plot_domain

    animation = @animate for frame in 1:nframes
        idx = max(2, round(Int, 1 + (length(xs) - 1) * (frame - 1) / (nframes - 1)))
        first_idx = max(1, idx - tail)

        plt = base_plot(tg.mesh)
        full_path_alpha = with_mesh ? 0.14 : 0.45
        full_path_width = with_mesh ? 0.40 : 0.55
        plot!(plt, xs, ys; linecolor="#AAB3BB", linealpha=full_path_alpha, linewidth=full_path_width)
        plot!(
            plt,
            xs[first_idx:idx],
            ys[first_idx:idx];
            linecolor="#111827",
            linealpha=0.95,
            linewidth=2.1,
        )
        scatter!(
            plt,
            [xs[idx]],
            [ys[idx]];
            markercolor="#E76F51",
            markerstrokecolor=:white,
            markerstrokewidth=1.2,
            markersize=4.0,
        )
    end

    gif(animation, joinpath(OUTDIR, output); fps=24)
end

save_static_assets(tg)
save_cycle_animation(tg, "cyclic_track_no_mesh.gif"; with_mesh=false)
save_cycle_animation(tg, "cyclic_track_with_mesh.gif"; with_mesh=true)
