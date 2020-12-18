using RayTracing
using Gridap
using Test

jsonfile = joinpath(@__DIR__,"pincell.json")
model = DiscreteModelFromFile(jsonfile)

@testset "Basic tests" begin

    tg = TrackGenerator(model, 8, 0.02)
    trace!(tg)
    segmentize!(tg)

    @testset "Tracing" begin
        @test tg.n_total_tracks == 420
        @test tg.n_tracks_x == [31, 74, 74, 31]
        @test tg.n_tracks_y == [74, 31, 31, 74]
        @test tg.n_tracks == [105, 105, 105, 105]
    end

    @testset "Azimuthal quadradure" begin
        @test isequal(tg.azimuthal_quadrature.n_azim, 8)
        @test isequal(tg.azimuthal_quadrature.n_azim_2, 4)
        @test isequal(tg.azimuthal_quadrature.n_azim_4, 2)
        @test isapprox(tg.azimuthal_quadrature.δ, 0.02)
        @test all(isapprox.(tg.azimuthal_quadrature.δs, 0.01994243696980254))
        @test tg.azimuthal_quadrature.ϕs ≈ [0.39670866289121387, 1.1740876639036828, 1.9675049896861103, 2.7448839906985794]
    end

    # @testset boundary points
end

@testset "Broken case" begin
    tg = TrackGenerator(model, 8, 0.2)
    trace!(tg)
    @test_broken segmentize!(tg) == error("This is an unexpected case. Please, submit an issue.")
end