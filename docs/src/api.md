```@meta
CurrentModule = RayTracing
```

# API Reference

## Main Workflow

```@docs
TrackGenerator
trace!
segmentize!
```

## Boundary Conditions

```@docs
RayTracing.BoundaryType
BoundaryConditions
reflective_boundaries
vacuum_boundaries
periodic_boundaries
RayTracing.get_boundary_condition_at
```

## Geometry Types

```@docs
Track
Segment
Mesh
```

## Geometry Helpers

```@docs
RayTracing.Point2D
RayTracing.advance_step
RayTracing.distance
RayTracing.midpoint
RayTracing.angle
RayTracing.is_approx
```

## Azimuthal Quadrature

```@docs
RayTracing.AzimuthalQuadrature
RayTracing.n_azim_total
RayTracing.n_azim_half
RayTracing.n_azim_quad
RayTracing.azimuthal_quadrant_1
RayTracing.azimuthal_quadrant_2
RayTracing.azimuthal_half_plane
RayTracing.supplementary_azimuthal_idx
RayTracing.init_weights!
```

## Solver Integration Internals

These helpers are useful for packages that consume RayTracing's track graph directly.

```@docs
RayTracing.bc_fwd
RayTracing.bc_bwd
RayTracing.dir_next_track_fwd
RayTracing.dir_next_track_bwd
RayTracing.universal_id
```
