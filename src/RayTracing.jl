# module CyclicRayTracing
module RayTracing

struct Segment
    xi
    xo
    ℓ

    τ

    element
end

struct Track
    id

    azim_idx

    ABC
    ϕ

    # alguna struct q sea un Point? StaticArrays.jl o Tensors.jl o GeometricalPredicates.jl
    xi
    xo
    ℓ

    segments::Vector{Segment}

    boundary_i
    boundary_o

    next_track_fwd
    next_track_bwd

    dir_next_track_fwd
    dir_next_track_bwd
end

# mmm la podriamos llamar a esta CyclicRayTracing o RayTracing tambien
struct Tracks
    mesh

    n_azim    # el numero de angulos azimutales en (0, 2π)
    n_azim_2  # el numero de angulos azimutales en (0, π)
    n_azim_4  # el numero de angulos azimutales en (0, π/2)

    dens
    spacing

    tiny_step

    n_tracks_x
    n_tracks_y
    n_tracks
    n_total_tracks

    track         # depende del indice del angulo azimutal y del numero de track
    track_by_id

    correct_volumes
    volumes

    τmax

    # quadrature (?)
end

# mesh es de gmsh, right? aun sin leer con GmshDiscreteModel
function Tracks(mesh, n_azim, spacing ; tiny_step=1e-8, trace=true)

    @assert n_azim > 0 "the number of azimuthal angles must be > 0."
    @assert n_azim % 4 == 0 "the number of azimuthal angles must be a multiple of 4."
    @assert spacing > 0 "the track spacing must be > 0."

    n_azim_2 = n_azim / 2
    n_azim_4 = n_azim / 4
    dens = inv(spacing)

    n_tracks_x = Vector{Int}(undef, n_azim_2)
    n_tracks_y = Vector{Int}(undef, n_azim_2)
    n_tracks   = Vector{Int}(undef, n_azim_2)

    # if !trace
    #     Tracks(mesh, n_azim, n_azim_2, n_azim_4, dens, spacing, tiny_step, n_tracks_x, n_tracks_y, n_tracks, 0, )
    # end

    # deltas entre los tracks en cada eje y perpendicular a los mismos
    T = typeof(dens)
    dx_eff = Vector{T}(undef, n_azim_2)
    dy_eff = Vector{T}(undef, n_azim_2)
    d_eff  = Vector{T}(undef, n_azim_2)

    # la bounding box la puedo buscar usando gmsh.model.getBoundingBox recorriendo todas las
    # entities y viendo cual tiene la bounding box mas grande, sin embargo aun no se que
    # objeto es el que estoy pasando aca como mesh, si el de GmshDiscreteModel o el de gmsh
    for entity in entities
        # mmm esto no es correcto si no tengo una entity que engloble a todas y tengo la
        # union de entities... me parece que no me queda otra que recorrer los nodos y
        # buscar el xmin, ymin, zmin, xmax, ymax, zmax y ya, como wasora:
        # https://github.com/seamplex/wasora/blob/acbbd98a08b62d8d3e1c709f8b48789d105f3589/src/mesh/mesh.c#L111
        # en este caso, creo que deberia usar ya el objecto GmshDiscreteModel, quizas, no lo
        # se
        current_bb = BoundingBox(gmsh.model.getBoundingBox(2, entity)...)
        bb = current_bb > bb ? current_bb : bb # hacer metodo > para bounding boxes?
    end

    width_x = bb.xmax - bb.xmin # dx(bb)
    width_y = bb.ymax - bb.ymin

    # barremos el primer cuadrante
    for i in 1:n_azim_4

        φ = π / n_azim_2 * (i - 1 / 2)

        n_tracks_x[i] = floor(width_x / spacing * abs(sin(φ))) + 1
        n_tracks_y[i] = floor(width_y / spacing * abs(cos(φ))) + 1
        n_tracks[i] = n_tracks_x[i] + n_tracks_y[i];

        # angulo efectivo, cerca del que deseamos
        ϕ = atan((width_y * n_tracks_x[i]) / (width_x * n_tracks_y[i]))
        # TODO: tambien hay q calcular cosas que iran luego a inicializar una quadrature
        # tracks->quadrature->azimuthal->phi[i] = ϕ

        # deltas efectivos
        dx_eff[i] = width_x / n_tracks_x[i]
        dy_eff[i] = width_y / n_tracks_y[i]
        d_eff[i] = dx_eff[i] * sin(ϕ)

        # idem arriba
        # tracks->quadrature->effective_spacing[i] = d_eff[i]

        # y los angulos suplementarios
        n_tracks_x[n_azim_2-i+1] = n_tracks_x[i]
        n_tracks_y[n_azim_2-i+1] = n_tracks_y[i]
        n_tracks[n_azim_2-i+1] = n_tracks[i]

        # idem arriba
        # tracks->quadrature->azimuthal->phi[tracks->n_azim_2-i+1] = π - ϕ

        dx_eff[n_azim_2-i+1] = dx_eff[i]
        dy_eff[n_azim_2-i+1] = dy_eff[i]
        d_eff[n_azim_2-i+1] = d_eff[i]

        # idem arriba
        # tracks->quadrature->effective_spacing[tracks->n_azim_2-i+1] = d_eff[i]
    end

    #  empiezo a allocar la "matriz" de tracks
    track = Vector{Vector{Track}}(undef, n_azim_2)

    # puntos de inicio y final de cada track (referenciados al (0,0) ficticio en el borde inferior izquierdo)
    # barro los dos primeros cuadrantes
    for i in 1:n_azim_2

        # recupero el angulo
        # ϕ = tracks->quadrature->azimuthal->phi[i]

        # alloco segun la cantidad de tracks dado el angulo i
        track[i] = Vector{Track}(undef, n_tracks[i])

        # puntos de inicio para cada track sobre el eje x
        for j in 1:n_tracks_x[i]

            # no podes referenciar a algo que aun no esta definido rami.
            # t = track[i][j]
            # entonces primero defino todas las cosas de un track y luego defino el track y lo
            # igualo a lo que corresponda

            # si el track apunta a la derecha o a la izquierda
            xi = i <= n_azim_4 ?
                Point(dx_eff[i] * (n_tracks_x[i] - j + 1 / 2), 0, 0) :
                Point(dx_eff[i] * (j - 1 / 2), 0, 0)




            track[i][j] = Track(...)
        end



        # al final pusheamos:
        track[i] = track_i
    end


    # return Tracks(mesh, n_azim, n_azim_2, n_azim_4, dens, spacing, tiny_step, n_tracks_x,
    #               n_tracks_y, n_tracks, )
end

struct BoundingBox{T}
    xmin::T
    ymin::T
    zmin::T
    xmax::T
    ymax::T
    zmax::T
end

BoundingBox(xmin, ymin, zmin, xmax, ymax, zmax) =
    BoundingBox(promote(xmin, ymin, zmin, xmax, ymax, zmax)...)

# hacer otros? ver bien syntax
import Base: >
function >(bb1::BoundingBox, bb2::BoundingBox)
    # xmin1, ymin1, zmin1, xmax1, ymax1, zmax1 = bb1
    # xmin2, ymin2, zmin2, xmax2, ymax2, zmax2 = bb2
end


# tambien se puede llamar `track`
function trace(tracks::Tracks) end


end
