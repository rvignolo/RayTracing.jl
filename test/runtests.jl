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

@testset "Reflection tests" begin

    bcs = BoundaryConditions(top=Vaccum, bottom=Reflective, left=Reflective, right=Vaccum)
    tg = TrackGenerator(model, 4, 0.8; bcs=bcs)
    trace!(tg)

    @testset "Track 1" begin
        track = tg.tracks_by_uid[1]
        next_track_fwd = track.next_track_fwd
        next_track_bwd = track.next_track_bwd
        DirNextTrackFwd = RayTracing.dir_next_track_fwd(track)
        DirNextTrackBwd = RayTracing.dir_next_track_bwd(track)
        @test RayTracing.bc_fwd(track) == Vaccum
        @test RayTracing.bc_bwd(track) == Reflective
        @test next_track_fwd.uid == 7
        @test next_track_bwd.uid == 6
        @test DirNextTrackFwd == RayTracing.Forward
        @test DirNextTrackBwd == RayTracing.Forward
    end

    @testset "Track 2" begin
        track = tg.tracks_by_uid[2]
        next_track_fwd = track.next_track_fwd
        next_track_bwd = track.next_track_bwd
        DirNextTrackFwd = RayTracing.dir_next_track_fwd(track)
        DirNextTrackBwd = RayTracing.dir_next_track_bwd(track)
        @test RayTracing.bc_fwd(track) == Vaccum
        @test RayTracing.bc_bwd(track) == Reflective
        @test next_track_fwd.uid == 8
        @test next_track_bwd.uid == 5
        @test DirNextTrackFwd == RayTracing.Forward
        @test DirNextTrackBwd == RayTracing.Forward
    end

    @testset "Track 3" begin
        track = tg.tracks_by_uid[3]
        next_track_fwd = track.next_track_fwd
        next_track_bwd = track.next_track_bwd
        DirNextTrackFwd = RayTracing.dir_next_track_fwd(track)
        DirNextTrackBwd = RayTracing.dir_next_track_bwd(track)
        @test RayTracing.bc_fwd(track) == Vaccum
        @test RayTracing.bc_bwd(track) == Reflective
        @test next_track_fwd.uid == 8
        @test next_track_bwd.uid == 5
        @test DirNextTrackFwd == RayTracing.Backward
        @test DirNextTrackBwd == RayTracing.Backward
    end

    @testset "Track 4" begin
        track = tg.tracks_by_uid[4]
        next_track_fwd = track.next_track_fwd
        next_track_bwd = track.next_track_bwd
        DirNextTrackFwd = RayTracing.dir_next_track_fwd(track)
        DirNextTrackBwd = RayTracing.dir_next_track_bwd(track)
        @test RayTracing.bc_fwd(track) == Vaccum
        @test RayTracing.bc_bwd(track) == Reflective
        @test next_track_fwd.uid == 7
        @test next_track_bwd.uid == 6
        @test DirNextTrackFwd == RayTracing.Backward
        @test DirNextTrackBwd == RayTracing.Backward
    end

    @testset "Track 5" begin
        track = tg.tracks_by_uid[5]
        next_track_fwd = track.next_track_fwd
        next_track_bwd = track.next_track_bwd
        DirNextTrackFwd = RayTracing.dir_next_track_fwd(track)
        DirNextTrackBwd = RayTracing.dir_next_track_bwd(track)
        @test RayTracing.bc_fwd(track) == Reflective
        @test RayTracing.bc_bwd(track) == Reflective
        @test next_track_fwd.uid == 3
        @test next_track_bwd.uid == 2
        @test DirNextTrackFwd == RayTracing.Forward
        @test DirNextTrackBwd == RayTracing.Forward
    end

    @testset "Track 6" begin
        track = tg.tracks_by_uid[6]
        next_track_fwd = track.next_track_fwd
        next_track_bwd = track.next_track_bwd
        DirNextTrackFwd = RayTracing.dir_next_track_fwd(track)
        DirNextTrackBwd = RayTracing.dir_next_track_bwd(track)
        @test RayTracing.bc_fwd(track) == Reflective
        @test RayTracing.bc_bwd(track) == Reflective
        @test next_track_fwd.uid == 4
        @test next_track_bwd.uid == 1
        @test DirNextTrackFwd == RayTracing.Forward
        @test DirNextTrackBwd == RayTracing.Forward
    end

    @testset "Track 7" begin
        track = tg.tracks_by_uid[7]
        next_track_fwd = track.next_track_fwd
        next_track_bwd = track.next_track_bwd
        DirNextTrackFwd = RayTracing.dir_next_track_fwd(track)
        DirNextTrackBwd = RayTracing.dir_next_track_bwd(track)
        @test RayTracing.bc_fwd(track) == Vaccum
        @test RayTracing.bc_bwd(track) == Vaccum
        @test next_track_fwd.uid == 4
        @test next_track_bwd.uid == 1
        @test DirNextTrackFwd == RayTracing.Backward
        @test DirNextTrackBwd == RayTracing.Backward
    end

    @testset "Track 8" begin
        track = tg.tracks_by_uid[8]
        next_track_fwd = track.next_track_fwd
        next_track_bwd = track.next_track_bwd
        DirNextTrackFwd = RayTracing.dir_next_track_fwd(track)
        DirNextTrackBwd = RayTracing.dir_next_track_bwd(track)
        @test RayTracing.bc_fwd(track) == Vaccum
        @test RayTracing.bc_bwd(track) == Vaccum
        @test next_track_fwd.uid == 3
        @test next_track_bwd.uid == 2
        @test DirNextTrackFwd == RayTracing.Backward
        @test DirNextTrackBwd == RayTracing.Backward
    end

    # a more complicated case
    tg = TrackGenerator(model, 8, 0.8; bcs=bcs)
    trace!(tg)

    @testset "Track 1" begin
        track = tg.tracks_by_uid[1]
        next_track_fwd = track.next_track_fwd
        next_track_bwd = track.next_track_bwd
        DirNextTrackFwd = RayTracing.dir_next_track_fwd(track)
        DirNextTrackBwd = RayTracing.dir_next_track_bwd(track)
        @test RayTracing.bc_fwd(track) == Vaccum
        @test RayTracing.bc_bwd(track) == Reflective
        @test next_track_fwd.uid == 11
        @test next_track_bwd.uid == 10
        @test DirNextTrackFwd == RayTracing.Forward
        @test DirNextTrackBwd == RayTracing.Forward
    end

    @testset "Track 2" begin
        track = tg.tracks_by_uid[2]
        next_track_fwd = track.next_track_fwd
        next_track_bwd = track.next_track_bwd
        DirNextTrackFwd = RayTracing.dir_next_track_fwd(track)
        DirNextTrackBwd = RayTracing.dir_next_track_bwd(track)
        @test RayTracing.bc_fwd(track) == Vaccum
        @test RayTracing.bc_bwd(track) == Reflective
        @test next_track_fwd.uid == 12
        @test next_track_bwd.uid == 10
        @test DirNextTrackFwd == RayTracing.Forward
        @test DirNextTrackBwd == RayTracing.Backward
    end

    @testset "Track 3" begin
        track = tg.tracks_by_uid[3]
        next_track_fwd = track.next_track_fwd
        next_track_bwd = track.next_track_bwd
        DirNextTrackFwd = RayTracing.dir_next_track_fwd(track)
        DirNextTrackBwd = RayTracing.dir_next_track_bwd(track)
        @test RayTracing.bc_fwd(track) == Vaccum
        @test RayTracing.bc_bwd(track) == Reflective
        @test next_track_fwd.uid == 12
        @test next_track_bwd.uid == 11
        @test DirNextTrackFwd == RayTracing.Backward
        @test DirNextTrackBwd == RayTracing.Backward
    end

    @testset "Track 4" begin
        track = tg.tracks_by_uid[4]
        next_track_fwd = track.next_track_fwd
        next_track_bwd = track.next_track_bwd
        DirNextTrackFwd = RayTracing.dir_next_track_fwd(track)
        DirNextTrackBwd = RayTracing.dir_next_track_bwd(track)
        @test RayTracing.bc_fwd(track) == Vaccum
        @test RayTracing.bc_bwd(track) == Reflective
        @test next_track_fwd.uid == 9
        @test next_track_bwd.uid == 8
        @test DirNextTrackFwd == RayTracing.Forward
        @test DirNextTrackBwd == RayTracing.Forward
    end

    @testset "Track 5" begin
        track = tg.tracks_by_uid[5]
        next_track_fwd = track.next_track_fwd
        next_track_bwd = track.next_track_bwd
        DirNextTrackFwd = RayTracing.dir_next_track_fwd(track)
        DirNextTrackBwd = RayTracing.dir_next_track_bwd(track)
        @test RayTracing.bc_fwd(track) == Vaccum
        @test RayTracing.bc_bwd(track) == Reflective
        @test next_track_fwd.uid == 9
        @test next_track_bwd.uid == 7
        @test DirNextTrackFwd == RayTracing.Backward
        @test DirNextTrackBwd == RayTracing.Forward
    end

    @testset "Track 6" begin
        track = tg.tracks_by_uid[6]
        next_track_fwd = track.next_track_fwd
        next_track_bwd = track.next_track_bwd
        DirNextTrackFwd = RayTracing.dir_next_track_fwd(track)
        DirNextTrackBwd = RayTracing.dir_next_track_bwd(track)
        @test RayTracing.bc_fwd(track) == Vaccum
        @test RayTracing.bc_bwd(track) == Reflective
        @test next_track_fwd.uid == 8
        @test next_track_bwd.uid == 7
        @test DirNextTrackFwd == RayTracing.Backward
        @test DirNextTrackBwd == RayTracing.Backward
    end

    @testset "Track 7" begin
        track = tg.tracks_by_uid[7]
        next_track_fwd = track.next_track_fwd
        next_track_bwd = track.next_track_bwd
        DirNextTrackFwd = RayTracing.dir_next_track_fwd(track)
        DirNextTrackBwd = RayTracing.dir_next_track_bwd(track)
        @test RayTracing.bc_fwd(track) == Reflective
        @test RayTracing.bc_bwd(track) == Reflective
        @test next_track_fwd.uid == 6
        @test next_track_bwd.uid == 5
        @test DirNextTrackFwd == RayTracing.Forward
        @test DirNextTrackBwd == RayTracing.Forward
    end

    @testset "Track 8" begin
        track = tg.tracks_by_uid[8]
        next_track_fwd = track.next_track_fwd
        next_track_bwd = track.next_track_bwd
        DirNextTrackFwd = RayTracing.dir_next_track_fwd(track)
        DirNextTrackBwd = RayTracing.dir_next_track_bwd(track)
        @test RayTracing.bc_fwd(track) == Vaccum
        @test RayTracing.bc_bwd(track) == Reflective
        @test next_track_fwd.uid == 6
        @test next_track_bwd.uid == 4
        @test DirNextTrackFwd == RayTracing.Backward
        @test DirNextTrackBwd == RayTracing.Forward
    end

    @testset "Track 9" begin
        track = tg.tracks_by_uid[9]
        next_track_fwd = track.next_track_fwd
        next_track_bwd = track.next_track_bwd
        DirNextTrackFwd = RayTracing.dir_next_track_fwd(track)
        DirNextTrackBwd = RayTracing.dir_next_track_bwd(track)
        @test RayTracing.bc_fwd(track) == Vaccum
        @test RayTracing.bc_bwd(track) == Vaccum
        @test next_track_fwd.uid == 5
        @test next_track_bwd.uid == 4
        @test DirNextTrackFwd == RayTracing.Backward
        @test DirNextTrackBwd == RayTracing.Backward
    end

    @testset "Track 10" begin
        track = tg.tracks_by_uid[10]
        next_track_fwd = track.next_track_fwd
        next_track_bwd = track.next_track_bwd
        DirNextTrackFwd = RayTracing.dir_next_track_fwd(track)
        DirNextTrackBwd = RayTracing.dir_next_track_bwd(track)
        @test RayTracing.bc_fwd(track) == Reflective
        @test RayTracing.bc_bwd(track) == Reflective
        @test next_track_fwd.uid == 2
        @test next_track_bwd.uid == 1
        @test DirNextTrackFwd == RayTracing.Forward
        @test DirNextTrackBwd == RayTracing.Forward
    end

    @testset "Track 11" begin
        track = tg.tracks_by_uid[11]
        next_track_fwd = track.next_track_fwd
        next_track_bwd = track.next_track_bwd
        DirNextTrackFwd = RayTracing.dir_next_track_fwd(track)
        DirNextTrackBwd = RayTracing.dir_next_track_bwd(track)
        @test RayTracing.bc_fwd(track) == Reflective
        @test RayTracing.bc_bwd(track) == Vaccum
        @test next_track_fwd.uid == 3
        @test next_track_bwd.uid == 1
        @test DirNextTrackFwd == RayTracing.Forward
        @test DirNextTrackBwd == RayTracing.Backward
    end

    @testset "Track 12" begin
        track = tg.tracks_by_uid[12]
        next_track_fwd = track.next_track_fwd
        next_track_bwd = track.next_track_bwd
        DirNextTrackFwd = RayTracing.dir_next_track_fwd(track)
        DirNextTrackBwd = RayTracing.dir_next_track_bwd(track)
        @test RayTracing.bc_fwd(track) == Vaccum
        @test RayTracing.bc_bwd(track) == Vaccum
        @test next_track_fwd.uid == 3
        @test next_track_bwd.uid == 2
        @test DirNextTrackFwd == RayTracing.Backward
        @test DirNextTrackBwd == RayTracing.Backward
    end
end
