using BenchmarkTools
using Gridap
using RayTracing

const SUITE = BenchmarkGroup()

const PINCELL_JSON = normpath(joinpath(@__DIR__, "..", "demo", "pincell.json"))
const MODEL = DiscreteModelFromFile(PINCELL_JSON)
const MESH = RayTracing.Mesh(MODEL)

function make_generator(n_azim, spacing, bcs, volume_correction)
    return TrackGenerator(MODEL, n_azim, spacing; bcs=bcs, volume_correction=volume_correction)
end

function traced_generator(n_azim, spacing, bcs, volume_correction)
    tg = make_generator(n_azim, spacing, bcs, volume_correction)
    trace!(tg)
    return tg
end

function segmented_generator(n_azim, spacing, bcs, volume_correction)
    tg = traced_generator(n_azim, spacing, bcs, volume_correction)
    segmentize!(tg)
    return tg
end

function first_track_segment(tg)
    for track in tg.tracks_by_uid
        isempty(track.segments) && continue
        return track, first(track.segments)
    end
    error("No segmented tracks are available in the benchmark model.")
end

function point_batch(n)
    return [RayTracing.Point2D(rand(), rand()) for _ in 1:n]
end

function shifted_point_batch(points)
    return [p + RayTracing.Point2D(rand() + 1, rand() + 1) for p in points]
end

function line_batch(points_a, points_b)
    return [
        RayTracing.general_form(points_a[i], points_b[i])
        for i in eachindex(points_a, points_b)
    ]
end

function advance_step_batch(points, steps, ϕs)
    s = 0.0
    @inbounds for i in eachindex(points, steps, ϕs)
        p = RayTracing.advance_step(points[i], steps[i], ϕs[i])
        s += p[1] + p[2]
    end
    return s
end

function distance_batch(points_a, points_b)
    s = 0.0
    @inbounds for i in eachindex(points_a, points_b)
        s += RayTracing.distance(points_a[i], points_b[i])
    end
    return s
end

function general_form_batch(points_a, points_b)
    s = 0.0
    @inbounds for i in eachindex(points_a, points_b)
        ABC = RayTracing.general_form(points_a[i], points_b[i])
        s += ABC[1] + ABC[2] + ABC[3]
    end
    return s
end

function element_volume_batch(mesh, cell_ids)
    s = 0.0
    @inbounds for cell_id in cell_ids
        s += RayTracing.element_volume(mesh, cell_id)
    end
    return s
end

function line_intersection_batch(lines_a, lines_b)
    s = 0.0
    @inbounds for i in eachindex(lines_a, lines_b)
        parallel, point = RayTracing.intersection(lines_a[i], lines_b[i])
        s += point[1] + point[2] + parallel
    end
    return s
end

const SAMPLE_TG = segmented_generator(8, 0.08, vacuum_boundaries(), false)
const SAMPLE_TRACK_SEGMENT = first_track_segment(SAMPLE_TG)
const SAMPLE_TRACK = SAMPLE_TRACK_SEGMENT[1]
const SAMPLE_SEGMENT = SAMPLE_TRACK_SEGMENT[2]
const SAMPLE_CELL = SAMPLE_SEGMENT.element
const SAMPLE_POINT = RayTracing.midpoint(SAMPLE_SEGMENT.p, SAMPLE_SEGMENT.q)
const BATCH_SIZE = 1024

const MESH_REF = Ref(MESH)
const SAMPLE_TRACK_REF = Ref(SAMPLE_TRACK)
const SAMPLE_POINT_REF = Ref(SAMPLE_POINT)

const GEOMETRY = SUITE["geometry"] = BenchmarkGroup()
const POINTS = GEOMETRY["points"] = BenchmarkGroup()
const MESHES = GEOMETRY["meshes"] = BenchmarkGroup()
const INTERSECTIONS = GEOMETRY["intersections"] = BenchmarkGroup()

POINTS["advance_step_batch"] = @benchmarkable advance_step_batch(points, steps, ϕs) setup=(
    points = point_batch($BATCH_SIZE);
    steps = rand($BATCH_SIZE);
    ϕs = rand($BATCH_SIZE)
)
POINTS["distance_batch"] = @benchmarkable distance_batch(points_a, points_b) setup=(
    points_a = point_batch($BATCH_SIZE);
    points_b = point_batch($BATCH_SIZE)
)
POINTS["general_form_batch"] = @benchmarkable general_form_batch(points_a, points_b) setup=(
    points_a = point_batch($BATCH_SIZE);
    points_b = shifted_point_batch(points_a)
)

MESHES["construct"] = @benchmarkable RayTracing.Mesh($MODEL) evals=1
MESHES["find_element"] = @benchmarkable RayTracing.find_element($(MESH_REF)[], $(SAMPLE_POINT_REF)[])
MESHES["element_volume_batch"] = @benchmarkable element_volume_batch(mesh, cell_ids) setup=(
    mesh = $(MESH_REF)[];
    cell_ids = rand(1:length(mesh.cell_nodes), $BATCH_SIZE)
)

INTERSECTIONS["line_line_batch"] = @benchmarkable line_intersection_batch(lines_a, lines_b) setup=(
    points_a = point_batch($BATCH_SIZE);
    points_b = shifted_point_batch(points_a);
    points_c = point_batch($BATCH_SIZE);
    points_d = shifted_point_batch(points_c);
    lines_a = line_batch(points_a, points_b);
    lines_b = line_batch(points_c, points_d)
)
INTERSECTIONS["track_cell"] = @benchmarkable RayTracing.intersections(
    $(MESH_REF)[], $SAMPLE_CELL, $(SAMPLE_TRACK_REF)[]
)

const WORKFLOWS = SUITE["workflows"] = BenchmarkGroup()
const CASES = (
    "coarse-vacuum" => (; n_azim=4, spacing=0.16, bcs=vacuum_boundaries(), volume_correction=false),
    "moderate-vacuum" => (; n_azim=8, spacing=0.08, bcs=vacuum_boundaries(), volume_correction=false),
    "moderate-reflective" => (; n_azim=8, spacing=0.08, bcs=reflective_boundaries(), volume_correction=false),
    "moderate-periodic" => (; n_azim=8, spacing=0.08, bcs=periodic_boundaries(), volume_correction=false),
    "volume-correction" => (; n_azim=8, spacing=0.08, bcs=vacuum_boundaries(), volume_correction=true),
)

for (name, case) in CASES
    group = WORKFLOWS[name] = BenchmarkGroup()
    n_azim = case.n_azim
    spacing = case.spacing
    bcs = case.bcs
    volume_correction = case.volume_correction

    group["construct"] = @benchmarkable make_generator(
        $n_azim, $spacing, $bcs, $volume_correction
    ) evals=1

    group["trace!"] = @benchmarkable trace!(tg) setup=(
        tg = make_generator($n_azim, $spacing, $bcs, $volume_correction)
    ) evals=1

    group["segmentize!"] = @benchmarkable segmentize!(tg) setup=(
        tg = traced_generator($n_azim, $spacing, $bcs, $volume_correction)
    ) evals=1

    group["segmentize! parallel"] = @benchmarkable segmentize!(tg; parallel=true) setup=(
        tg = traced_generator($n_azim, $spacing, $bcs, $volume_correction)
    ) evals=1

    group["trace+segmentize!"] = @benchmarkable begin
        trace!(tg)
        segmentize!(tg)
    end setup=(
        tg = make_generator($n_azim, $spacing, $bcs, $volume_correction)
    ) evals=1

    group["trace+segmentize! parallel"] = @benchmarkable begin
        trace!(tg)
        segmentize!(tg; parallel=true)
    end setup=(
        tg = make_generator($n_azim, $spacing, $bcs, $volume_correction)
    ) evals=1
end
