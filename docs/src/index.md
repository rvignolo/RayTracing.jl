```@meta
CurrentModule = RayTracing
```

# RayTracing.jl

RayTracing.jl generates two-dimensional Method of Characteristics tracks over Gridap
unstructured meshes. It is designed as the geometry and tracking layer used by transport
codes such as [NeutronTransport.jl](https://github.com/rvignolo/NeutronTransport.jl).

The package builds straight tracks for azimuthal quadrature directions, connects them
according to rectangular-domain boundary conditions, segments each track by mesh element,
and stores per-cell ray-tracing volumes for downstream transport sweeps.

## Workflow

1. Load or build a `Gridap.Geometry.UnstructuredDiscreteModel`.
2. Construct a [`TrackGenerator`](@ref) with an azimuthal angle count and target spacing.
3. Run [`trace!`](@ref) to create tracks and connect boundary crossings.
4. Run [`segmentize!`](@ref) to split tracks into cell-local [`Segment`](@ref)s.
5. Use `tg.tracks_by_uid`, `tg.volumes`, and the azimuthal quadrature data in a solver.

Track segmentation uses `parallel=:auto` by default. When Julia is running with multiple
threads, `segmentize!(tg)` uses threaded segmentation.

## Quick Start

```julia
using Gridap
using RayTracing

model = DiscreteModelFromFile("pincell.json")

n_azim = 8
spacing = 0.08
bcs = reflective_boundaries()

tg = TrackGenerator(model, n_azim, spacing; bcs)
trace!(tg)
segmentize!(tg)

length(tg.tracks_by_uid)
sum(length(track.segments) for track in tg.tracks_by_uid)
```

## Boundary Conditions

RayTracing supports rectangular-domain vacuum, reflective, and periodic boundaries:

```julia
vacuum_boundaries()
reflective_boundaries()
periodic_boundaries()

BoundaryConditions(
    top=Vacuum,
    bottom=Reflective,
    left=Periodic,
    right=Periodic,
)
```

## Plotting

RayTracing ships lightweight Plots.jl recipes through RecipesBase. Plots.jl is not required
to use the core package, but if it is available you can inspect the generated geometry:

```julia
using Plots

plot(tg.mesh)
plot(tg)
plot(tg.tracks_by_uid)
```

The README media can be regenerated with:

```bash
julia --project=demo demo/pincell.jl
```

## Index

```@index
```
