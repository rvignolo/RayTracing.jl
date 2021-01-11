using RayTracing
using Gridap
using Test

jsonfile = joinpath(@__DIR__,"../demo/pincell.json")
model = DiscreteModelFromFile(jsonfile)

@testset "Principal tests" begin

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
        @test isequal(RayTracing.nazim(tg.azimuthal_quadrature), 8)
        @test isequal(RayTracing.nazim2(tg.azimuthal_quadrature), 4)
        @test isequal(RayTracing.nazim4(tg.azimuthal_quadrature), 2)
        @test isapprox(tg.azimuthal_quadrature.δ, 0.02)
        @test all(isapprox.(tg.azimuthal_quadrature.δs, 0.01994243696980254))
        @test tg.azimuthal_quadrature.ϕs ≈ [0.39670866289121387, 1.1740876639036828, 1.9675049896861103, 2.7448839906985794]
    end

    @testset "Entry and exit points" begin
        for track in tg.tracks_by_uid
            @test isapprox(track.p, track.segments[begin].p)
            @test isapprox(track.q, track.segments[end].q)
        end
    end

    @testset "Track length" begin
        for track in tg.tracks_by_uid
            l1 = track.ℓ
            l2 = sum(RayTracing.ℓ.(track.segments))
            @test isapprox(l1, l2)
        end
    end
end
