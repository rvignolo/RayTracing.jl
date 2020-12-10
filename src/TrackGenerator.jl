
# mmm la podriamos llamar a esta CyclicRayTracing o RayTracing tambien
# openMoC la llama TrackGenerator
struct Tracks{M,Q<:Quadrature,T<:Real}
    mesh::M
    # bb::BoundingBox

    quadrature::Q

    n_tracks_x::Vector{Int}
    n_tracks_y::Vector{Int}
    n_tracks::Vector{Int}
    n_total_tracks::Int

    tracks::Vector{Vector{Track}}   # depende del indice del angulo azimutal y del numero de track
    tracks_by_gid::Vector{Track}

    tiny_step::T
    correct_volumes::Bool
    volumes::Vector{T}

    # aunque un metodo seria muy costoso, asi que por ahi si convenga storearlo a medida q calculamos cada segmento
    # τmax # me parece que esto se deberia calcular con un metodo, recorriendo los τ de los segments
end

# en realidad mesh se suele llamar model (UnstructuredDiscreteModel), a menos que mandemos
# el objeto de gmsh, que no creo
function Tracks(mesh, n_azim, δ::T; tiny_step::T=1e-8, correct_volumes=false) where {T<:Real}

    @assert n_azim > 0 "the number of azimuthal angles must be > 0."
    @assert iszero(rem(n_azim, 4)) "the number of azimuthal angles must be a multiple of 4."
    @assert δ > 0 "the azimuthal spacing must be > 0."

    n_azim_2 = div(n_azim, 2)
    n_azim_4 = div(n_azim, 4)

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

        n_tracks_x[i] = floor(width_x / δ * abs(sin(φ))) + 1
        n_tracks_y[i] = floor(width_y / δ * abs(cos(φ))) + 1
        n_tracks[i] = n_tracks_x[i] + n_tracks_y[i]

        # los suplementarios tienen la misma cantidad
        j = n_azim_2 - i + 1
        n_tracks_x[j] = n_tracks_x[i]
        n_tracks_y[j] = n_tracks_y[i]
        n_tracks[j] = n_tracks[i]
    end

    n_total_tracks = sum(n_tracks)

    tracks = Vector{Vector{Track}}(undef, n_azim_2)
    for i in 1:n_azim_2
        tracks[i] = Vector{Track}(undef, n_tracks[i])
    end

    tracks_by_gid = Vector{Track}(undef, n_total_tracks)

    # ncells = get_number_cells(mesh) # buscar como se hace
    ncells = 100
    volumes = Vector{T}(undef, ncells)

    δs = Vector{T}(undef, n_azim_2)
    ϕs = Vector{T}(undef, n_azim_2)
    ωₐ = Vector{T}(undef, n_azim_2)

    # por ahi tambien podriamos permitir que se defina aqui la polar quadrature, con inputs
    n_polar = n_polar_2 = 0
    sinθs = Vector{T}(undef, n_polar)
    θs = Vector{T}(undef, n_polar)
    ωₚ = Vector{T}(undef, n_polar)

    polarquad = PolarQuadrature(n_polar, n_polar_2, sinθs, θs, ωₚ)

    ω = Matrix{T}(undef, n_azim_2, n_polar_2)

    quadrature = Quadrature(n_azim, n_azim_2, n_azim_4, δ, δs, ϕs, ωₐ, polarquad, ω)

    return Tracks(
        mesh, quadrature, n_tracks_x, n_tracks_y, n_tracks, n_total_tracks, tracks,
        tracks_by_gid, tiny_step, correct_volumes, volumes
    )
end

function trace(t::Tracks)
    @unpack quadrature = t
    @unpack n_tracks_x, n_tracks_y, n_tracks = t
    @unpack tracks, tracks_by_gid = t # mmm...
    @unpack n_azim_2, n_azim_4, δs, ϕs = quadrature

    # effective azimuthal spacings, used for some computations (not stored)
    T = eltype(δs)
    δx = Vector{T}(undef, n_azim_2)
    δy = Vector{T}(undef, n_azim_2)

    # @unpack bb = tracks
    # width_x = bb.xmax - bb.xmin # dx(bb)
    # width_y = bb.ymax - bb.ymin
    width_x = 5.
    width_y = 5.

    # first quadrant
    for i in 1:n_azim_4

        # effective azimuthal angle
        ϕ = ϕs[i] = atan((width_y * n_tracks_x[i]) / (width_x * n_tracks_y[i]))

        # effective azimuthal spacings
        δx[i] = width_x / n_tracks_x[i]
        δy[i] = width_y / n_tracks_y[i]
        δs[i] = δx[i] * sin(ϕ)

        # suplementarios:
        j = n_azim_2 - i + 1
        ϕs[j] = π - ϕ
        δx[j] = δx[i]
        δy[j] = δy[i]
        δs[j] = δs[i]
    end

    for i in 1:n_azim_2

        # get azimuthal angle
        ϕ = ϕs[i]

        # empiezo con los tracks que se originan sobre el eje x, que pueden apuntar a la
        # derecha o izquierda y , que pueden salir por y = width_y o x = width_x o x = 0
        for j in 1:n_tracks_x[i]

            # calculo el origin (dependiendo si apunta a derecha o izquierda)
            # hacer una funcion que indique si un track apunta a derecha o izquierda, aunque
            # aca no la puedo usar porque aun no tengo definido el Track
            if i <= n_azim_4
                xi = Point(δx[i] * (n_tracks_x[i] - j + 1 / 2), 0)
            else
                xi = Point(δx[i] * (j - 1 / 2), 0)
            end

            # calculo el end
            m = tan(ϕ) # useful

            # independientemente si apunta a derecha o izquierda puede salir por y = width_y
            xo = Point(xi[1] - (xi[2] - width_y) / m, width_y)

            # si no es punto de salida, seguimos buscando
            if !(0 <= xo[1] <= width_x)

                # aquellos que apuntan a la derecha, pueden salir tambien por x = width_x
                if i <= n_azim_4
                    xo = Point(width_x, xi[2] + m * (width_x - xi[1]))
                # aquellos que apuntan a la izquierda, pueden salir tambien por x = 0
                else
                    xo = Point(0, xi[2] - m * xi[1])
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
                xi = Point(0, δy[i] * (j - 1 / 2))
            else
                xi = Point(width_x, δy[i] * (j - 1 / 2))
            end

            ####################################################
            # FIXME: SPOT rule! esto se repite con lo de arriba

            # calculo el end
            m = tan(ϕ) # useful
            # independientemente si apunta a derecha o izquierda puede salir por y = width_y
            xo = Point(xi[1] - (xi[2] - width_y) / m, width_y)

            # si no es punto de salida, seguimos buscando
            if !(0 <= xo[1] <= width_x)

                # aquellos que apuntan a la derecha, pueden salir tambien por x = width_x
                if i <= n_azim_4
                    xo = Point(width_x, xi[2] + m * (width_x - xi[1]))
                # aquellos que apuntan a la izquierda, pueden salir tambien por x = 0
                else
                    xo = Point(0, xi[2] - m * xi[1])
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
    end

    id = 1
    for i in 1:n_azim_2, j in 1:n_tracks[i]
        tracks_by_gid[id] = tracks[i][j]
        id += 1
    end


    # defino los segmentos:
    #   1. tienen un origin
    #   2. tienen un end
    #   3. tienen una longitud
    #   4. tienen asociado un elemento
    #   5. dado el elemento, tienen un tau que podemos computar

    return t
end