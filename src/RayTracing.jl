# module CyclicRayTracing
module RayTracing

# IDEAS:
# 1. GridapEmbedded.jl para describir geometrias como open moc, constructive solid geometry.
#    Aca podria terminar implementando como tiene open moc las intersections entre figuras y
#    tracks.
# 2. mirar Meshes.jl y a la organizacion que pertenece.
# 3. anotar todas las funciones que necesito, como mesh_find_element(mesh, x) de wasora
#    (encontrar la cell que tiene el point), etc, y consultar si existen.


struct Segment{I,O,T}
    xi::I # son Gridap.Point o podrian ser un SVector ya que es semejante y tienen banda de operaciones hechas
    xo::O
    ℓ::T
    τ::T
    # element::E
end

# ID es id y AZMI es azim_idx
# struct Track{ID,AZMI,I,O,L,T}
struct Track{T,I,O,S}
    ϕ::T
    xi::I
    xo::O
    ℓ::T
    n::S
    segments::Vector{Segment}

    # por ahora comento esto:

    # boundary_i::PhysicalEntity
    # boundary_o::PhysicalEntity

    # next_track_fwd::Track
    # next_track_bwd::Track

    # dir_next_track_fwd::Int
    # dir_next_track_bwd::Int
end

# mmm la podriamos llamar a esta CyclicRayTracing o RayTracing tambien
struct Tracks{M,T}
    mesh::M

    # bb::BoundingBox

    n_azim::Int    # el numero de angulos azimutales en (0, 2π)
    n_azim_2::Int  # el numero de angulos azimutales en (0, π)
    n_azim_4::Int  # el numero de angulos azimutales en (0, π/2)

    dens::T
    spacing::T
    tiny_step::T

    n_tracks_x::Vector{Int}
    n_tracks_y::Vector{Int}
    n_tracks::Vector{Int}
    n_total_tracks::Int

    tracks::Vector{Vector{Track}}   # depende del indice del angulo azimutal y del numero de track
    tracks_by_id::Vector{Track}

    correct_volumes::Bool
    volumes::Vector{T}

    # τmax # me parece que esto se deberia calcular con un metodo, recorriendo los τ de los segments

    # quadrature (?)
end

# en realidad mesh se suele llamar model (UnstructuredDiscreteModel), a menos que mandemos
# el objeto de gmsh, que no creo
function Tracks(mesh, n_azim, spacing; tiny_step=1e-8, correct_volumes=false)

    @assert n_azim > 0 "the number of azimuthal angles must be > 0."
    @assert iszero(rem(n_azim, 4)) "the number of azimuthal angles must be a multiple of 4."
    @assert spacing > 0 "the track spacing must be > 0."

    n_azim_2 = div(n_azim, 2)
    n_azim_4 = div(n_azim, 4)

    dens = inv(spacing)

    n_tracks_x = Vector{Int}(undef, n_azim_2)
    n_tracks_y = Vector{Int}(undef, n_azim_2)
    n_tracks   = Vector{Int}(undef, n_azim_2)

    # obtenener bounding box usando get_nodal_coordinates(model)
    # bb = bounding_box(mesh) # BoundingBox(mesh)

    # width_x = bb.xmax - bb.xmin # dx(bb)
    # width_y = bb.ymax - bb.ymin
    width_x = 5.
    width_y = 5.

    # recorro un cuadrante y armo los suplementarios tambien
    for i in 1:n_azim_4
        φ = π / n_azim_2 * (i - 1 / 2)

        n_tracks_x[i] = floor(width_x / spacing * abs(sin(φ))) + 1
        n_tracks_y[i] = floor(width_y / spacing * abs(cos(φ))) + 1
        n_tracks[i] = n_tracks_x[i] + n_tracks_y[i]

        # los suplementarios tienen la misma cantidad
        n_tracks_x[n_azim_2-i+1] = n_tracks_x[i]
        n_tracks_y[n_azim_2-i+1] = n_tracks_y[i]
        n_tracks[n_azim_2-i+1] = n_tracks[i]
    end

    n_total_tracks = sum(n_tracks)

    track = Vector{Vector{Track}}(undef, n_azim_2)
    for i in 1:n_azim_2
        # alloco segun la cantidad de tracks dado el angulo i
        track[i] = Vector{Track}(undef, n_tracks[i])
    end

    track_by_id = Vector{Track}(undef, n_total_tracks)
    # mas adelante cuando ya tengamos asignados cada Track, hare:
    # track_by_id[id] = track[i][j] # con id computado
    # y esto seguro ya involucra que es un pointer (si modifico uno, modifica el otro)

    # ncells = get_number_cells(mesh) # buscar como se hace
    ncells = 100
    volumes = Vector{typeof(dens)}(undef, ncells)

    return Tracks(
        mesh, n_azim, n_azim_2, n_azim_4, dens, spacing, tiny_step, n_tracks_x, n_tracks_y,
        n_tracks, n_total_tracks, track, track_by_id, correct_volumes, volumes
    )
end

function trace(t::Tracks)

    @unpack dens, n_azim_2, n_azim_4 = t
    @unpack n_tracks_x, n_tracks_y, n_tracks = t
    @unpack tracks, tracks_by_id = t # mmm...

    # auxiliary
    T = typeof(dens)
    dx_eff = Vector{T}(undef, n_azim_2)
    dy_eff = Vector{T}(undef, n_azim_2)
    d_eff  = Vector{T}(undef, n_azim_2)

    # @unpack bb = tracks
    # width_x = bb.xmax - bb.xmin # dx(bb)
    # width_y = bb.ymax - bb.ymin
    width_x = 5.
    width_y = 5.

    # inicializo quadrature, q aun no la tengo definida dentro de Tracks
    for i in 1:n_azim_4

        # angulo efectivo, cerca del que deseamos
        ϕ = atan((width_y * n_tracks_x[i]) / (width_x * n_tracks_y[i]))

        # TODO: tambien hay q calcular cosas que iran luego a inicializar una quadrature
        # tracks->quadrature->azimuthal->phi[i] = ϕ

        # deltas efectivos
        dx_eff[i] = width_x / n_tracks_x[i]
        dy_eff[i] = width_y / n_tracks_y[i]
        d_eff[i] = dx_eff[i] * sin(ϕ)

        # tracks->quadrature->effective_spacing[i] = d_eff[i]

        # suplementarios:
        # tracks->quadrature->azimuthal->phi[tracks->n_azim_2-i+1] = π - ϕ
        dx_eff[n_azim_2-i+1] = dx_eff[i]
        dy_eff[n_azim_2-i+1] = dy_eff[i]
        d_eff[n_azim_2-i+1] = d_eff[i]
        # tracks->quadrature->effective_spacing[tracks->n_azim_2-i+1] = d_eff[i]
    end

    for i in 1:n_azim_2

        # recupero el angulo
        # ϕ = tracks->quadrature->azimuthal->phi[i]
        ϕ = atan((width_y * n_tracks_x[i]) / (width_x * n_tracks_y[i]))
        ϕ = i > n_azim_4 ? π - ϕ : ϕ

        # empiezo con los tracks que se originan sobre el eje x, que pueden apuntar a la
        # derecha o izquierda y , que pueden salir por y = width_y o x = width_x o x = 0
        for j in 1:n_tracks_x[i]

            # calculo el origin (dependiendo si apunta a derecha o izquierda)
            # hacer una funcion que indique si un track apunta a derecha o izquierda, aunque
            # aca no la puedo usar porque aun no tengo definido el Track
            if i <= n_azim_4
                xi = Gridap.Point(dx_eff[i] * (n_tracks_x[i] - j + 1 / 2), 0)
            else
                xi = Gridap.Point(dx_eff[i] * (j - 1 / 2), 0)
            end

            # calculo el end
            m = tan(ϕ) # useful

            # independientemente si apunta a derecha o izquierda puede salir por y = width_y
            xo = Gridap.Point(xi[1] - (xi[2] - width_y) / m, width_y)

            # si no es punto de salida, seguimos buscando
            if !(0 <= xo[1] <= width_x)

                # aquellos que apuntan a la derecha, pueden salir tambien por x = width_x
                if i <= n_azim_4
                    xo = Gridap.Point(width_x, xi[2] + m * (width_x - xi[1]))
                # aquellos que apuntan a la izquierda, pueden salir tambien por x = 0
                else
                    xo = Gridap.Point(0, xi[2] - m * xi[1])
                end

                if !(0 <= xo[2] <= width_y)
                    error("track exit point was not found.")
                end
            end

            ℓ = norm(xi - xo)
            # ℓ = distance(xi, xo) # que libreria la tiene, Gridap quizas

            # TODO: mejor calcularlo con LinearAlgebra y aca falta normalizar! llamarlo n?
            A = xi[2] - xo[2]
            B = xo[1] - xi[1]
            C = xi[1] * xo[2] - xo[1] * xi[2] # esto pareciera ser el determinante
            ABC = SVector(A, B, C) # por ahi tambien podria ser un Point...
            n = ABC / norm(ABC)

            segments = Vector{Segment}(undef, 0)

            tracks[i][j] = Track(ϕ, xi, xo, ℓ, n, segments)
        end

        # recorremos aquellos que inician en el eje y
        for j in 1:n_tracks_y[i]

            # aquellos que apuntan a la derecha inician en x = 0 y los que apuntan a la
            # izquierda inician en x = width_x
            if i <= n_azim_4
                xi = Gridap.Point(0, dy_eff[i] * (j - 1 / 2))
            else
                xi = Gridap.Point(width_x, dy_eff[i] * (j - 1 / 2))
            end

            ####################################################
            # FIXME: SPOT rule! esto se repite con lo de arriba

            # calculo el end
            m = tan(ϕ) # useful
            # independientemente si apunta a derecha o izquierda puede salir por y = width_y
            xo = Gridap.Point(xi[1] - (xi[2] - width_y) / m, width_y)

            # si no es punto de salida, seguimos buscando
            if !(0 <= xo[1] <= width_x)

                # aquellos que apuntan a la derecha, pueden salir tambien por x = width_x
                if i <= n_azim_4
                    xo = Gridap.Point(width_x, xi[2] + m * (width_x - xi[1]))
                # aquellos que apuntan a la izquierda, pueden salir tambien por x = 0
                else
                    xo = Gridap.Point(0, xi[2] - m * xi[1])
                end

                if !(0 <= xo[2] <= width_y)
                    error("track exit point was not found.")
                end
            end

            ℓ = norm(xi - xo)
            # ℓ = distance(xi, xo) # que libreria la tiene, Gridap quizas

            # TODO: mejor calcularlo con LinearAlgebra y aca falta normalizar! llamarlo n?
            A = xi[2] - xo[2]
            B = xo[1] - xi[1]
            C = xi[1] * xo[2] - xo[1] * xi[2] # esto pareciera ser el determinante
            ABC = SVector(A, B, C) # por ahi tambien podria ser un Point...
            n = ABC / norm(ABC)

            segments = Vector{Segment}(undef, 0)
            ####################################################

            tracks[i][n_tracks_x[i]+j] = Track(ϕ, xi, xo, ℓ, n, segments)
        end

        id = 1
        for i in 1:n_azim_2, j in 1:n_tracks[i]
            tracks_by_id[id] = tracks[i][j]
            id += 1
        end


    end


    # defino los segmentos:
    #   1. tienen un origin
    #   2. tienen un end
    #   3. tienen una longitud
    #   4. tienen asociado un elemento
    #   5. dado el elemento, tienen un tau que podemos computar

end


