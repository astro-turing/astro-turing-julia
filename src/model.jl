# ~\~ language=Julia filename=src/model.jl
# ~\~ begin <<lit/lheureux-model.md|src/model.jl>>[init]
using Printf: @printf
using LinearAlgebra: Tridiagonal

struct Param
    μ_a :: Float64
    μ_w :: Float64
    ρ_a :: Float64
    ρ_c :: Float64
    ρ_t :: Float64
    ρ_w :: Float64
    d0_ca :: Float64
    d0_co3 :: Float64
    β :: Float64
    b :: Float64
    k_a :: Float64
    k_c :: Float64
    k_1 :: Float64
    k_2 :: Float64
    k_3 :: Float64
    k_4 :: Float64
    m :: Float64
    n :: Float64
    x_d :: Float64
    h_d :: Float64
    l :: Float64
    s :: Float64
    c0_a :: Float64
    c0_c :: Float64
    s0_ca :: Float64
    s0_co3 :: Float64
    ϕ_inf :: Float64           # set to 0.01, see p. 7 bottom left
    g :: Float64               # gravity, set to 9.81 m/s^2
    ϕ_in :: Float64            # initial porosity
    ϕ0 :: Float64
    ϵ :: Float64
end

struct State
    c_a :: Vector{Float64}
    c_c :: Vector{Float64}
    s_ca :: Vector{Float64}
    s_co3 :: Vector{Float64}
    ϕ :: Vector{Float64}
end

SCENARIO_A(ϕ_in :: Float64, ϕ0 :: Float64) = Param(
    100.09, 18.45, 2.950, 2.710, 2.8, 1.023, 131.9, 272.6, 0.1, 5, 10^(-6.19), 10^(-6.37),
    1.0, 1.0, 0.1, 0.1, 2.48, 2.80, 50, 100, 500, 0.1, 0.6, 0.3, 0.326e-3, 0.326e-3, 0.01, 9.81,
    ϕ_in, ϕ0, 0.01)

SCENARIO_B(ϕ_in :: Float64, ϕ0 :: Float64) = Param(
    100.09, 18.45, 2.950, 2.710, 2.8, 1.023, 131.9, 272.6, 0.01, 10, 10^(-6.19), 10^(-6.37),
    0.01, 0.01, 0.001, 0.001, 2.48, 2.80, 50, 100, 500, 0.01, 0.6, 0.3, 0.326e-3, 0.326e-3, 0.01, 9.81,
    ϕ_in, ϕ0, 0.01)

function diffusion_matrix(c_d :: Vector{T}) where T <: Real
    Tridiagonal(c_d[2:end], -2 .* c_d, c_d[1:end-1])
end

function advection_matrix(c_a :: Union{T,Vector{T}}, sigma :: Vector{T}) where T <: Real
    Tridiagonal((c_a.*(sigma .+ 1))[2:end],
                2 .* c_a .* sigma,
                (c_a.*(sigma .- 1))[1:end-1])
end

"""
Finite difference of an array `y` in an upwind scheme given velocities `a`
"""
function upwind_dy(y::Vector{T}, a::Vector{T}) where T <: Real
    dy = Array{T}(undef, length(y))
    for (i, a) in enumerate(a)
            if i == 1 || i == length(y)
                    dy[i] = 0
            else
                    dy[i] = a < 0 ? a * (y[i+1] - y[i]) : a * (y[i] - y[i-1])
            end
    end
    dy
end

function clip(x, a, b)
    x < a ? a : (x > b ? b : x)
end

function clip_to_pos(x)
    x < 0 ? 0 : x
end

function clip_ara_cal(ara, cal)
    for (a, c) in zip(ara, cal)
        if (1 - a - c < 0)
            a = 1 - c
        end
        a = clip(a, 0, 1)
        c = clip(c, 0, 1)
    end
end

function propagators(p::Param, Δt::Float64, Δx::Float64)
    da = p.k_2 * p.d0_ca / p.s^2
    λ = p.k_3 / p.k_2
    ν_1 = p.k_1 / p.k_2
    ν_2 = p.k_4 / p.k_3

    # Equations 15: K
    conductivity(ϕ::Float64) = p.β * ϕ^3 / (1 - ϕ)^2 * conductivity_correction(ϕ)

    # Equation 17: F
    conductivity_correction(ϕ::Float64) = 1 - exp(10 - 10/ϕ)

    # Equation 24: H
    elastic_response(ϕ_nr::Float64) = -1 / (p.b * (ϕ_nr - p.ϕ_inf))

    # Equation 25: D_ϕ, this function is taken as a constant in the diffusion terms
    # and is completely ignored in the advection terms.
    porosity_diffusion(ϕ::Float64) =
        (p.β / p.b * p.g * p.ρ_w) * ϕ^3 / (1 - ϕ) / (p.ϕ_in - p.ϕ_inf) *
        conductivity_correction(ϕ)

    # This is strictly speaking a function of ϕ, but taken to be constant
    d_ϕ = porosity_diffusion(p.ϕ_in) / p.d0_ca
    x_star = p.d0_ca / p.s
    t_star = p.d0_ca / p.s^2

    function f_θ(x::Float64)
        x_prime = x * x_star
        if p.x_d < x_prime && x_prime < (p.x_d + p.h_d)
            1
        else
            0
        end
    end

    x = collect(LinRange(0.0, 1.0, 1024))
    θ = f_θ.(x)

    # § 2.5 states that ρ_s is taken as a constant, therefore delta will be a constant
    ρ_s = p.ρ_t / (1 - p.c0_a * (1 - p.ρ_t / p.ρ_a) - p.c0_c * (1 - p.ρ_t / p.ρ_c))
    δ = ρ_s / (p.μ_a * sqrt(p.k_c))

    # Functions of state
    #######################################################
    d_ca(s::State) = 1.0 ./ (1 .- 2 .* log.(s.ϕ))
    d_co3(s::State) = (p.d0_co3 / p.d0_ca) ./ (1 .- 2 .* log.(s.ϕ))

    sqr(x) = x*x
    # Equations 46 & 47: TODO reduce common subexpressions
    velocity_u(s::State) =
        1 - conductivity(p.ϕ0) / p.s * (1 - p.ϕ0) * (ρ_s / p.ρ_w - 1) .+
        conductivity.(s.ϕ) ./ p.s .* (1 .- s.ϕ) .* (ρ_s / p.ρ_w - 1)
    velocity_w(s::State) =
        1 - conductivity(p.ϕ0) / p.s * (1 - p.ϕ0) * (ρ_s / p.ρ_w - 1) .+
        conductivity.(s.ϕ) ./ p.s .* sqr.(1 .- s.ϕ) ./ s.ϕ .* (ρ_s / p.ρ_w - 1)

    Ω_pa(s::State) = clip_to_pos.(s.s_ca .* s.s_co3 .* (p.k_c / p.k_a) .- 1).^(p.m)
    Ω_da(s::State) = clip_to_pos.(1 .- s.s_ca .* s.s_co3 .* (p.k_c / p.k_a)).^(p.m) .* θ
    Ω_pc(s::State) = clip_to_pos.(s.s_ca .* s.s_co3 .- 1).^(p.n)
    Ω_dc(s::State) = clip_to_pos.(1 .- s.s_ca .* s.s_co3).^(p.n)

    ξ_1(s::State) =
        (1 .- s.c_a) .* s.c_a .* (Ω_da(s) .- ν_1 .* Ω_pa(s)) .+
        λ .* s.c_a .* s.c_c .* (Ω_pc(s) - ν_2 .* Ω_dc(s))
    ξ_2(s::State) = 
        λ .* (1 .- s.c_c) .* s.c_c .* (Ω_pc(s) - ν_2 .* Ω_dc(s)) .+
        s.c_a .* s.c_c .* (Ω_da(s) - ν_1 .* Ω_pa(s))
    ξ_3(s::State) =
        s.c_a .* (Ω_da(s) .- ν_1 .* Ω_pa(s)) .-
        λ .* s.c_c .* (Ω_pc(s) .- ν_2 .* Ω_dc(s))

    function sigma_peclet(w, d)
        c_f = (Δx/2) .* w ./ d
        coth.(c_f) - 1 ./ c_f
    end

    # PDE
    ###########################################################
    function half_step(s::State)
        u = velocity_u(s)
        w = velocity_w(s)

        c_a_half = s.c_a + (Δt/(2*Δx)) .* upwind_dy(s.c_a, u) - (da * Δt/2) .* ξ_1(s)
        c_c_half = s.c_a + (Δt/(2*Δx)) .* upwind_dy(s.c_c, u) - (da * Δt/2) .* ξ_2(s)
        clip_ara_cal(c_a_half, c_c_half)

        ϕ_half = clamp.(
            s.ϕ - (Δt/(2*Δx))   .* (advection_matrix(-1.0, sigma_peclet(w, d_ϕ)) * (s.ϕ .* w)) .+ 
                  (Δt/(4*Δx^2)) .* diffusion_matrix(fill(d_ϕ, 1024)) .* s.ϕ .+
                  (Δt*da/2) .* (1 .- s.ϕ) .* ξ_3(s),
            p.ϵ, 1 - p.ϵ)
        s_ca_half = ifelse.(
            s.ϕ .<= p.ϵ,
            s.s_ca - (Δt/(2*Δx)) .* (advection_matrix(w, sigma_peclet(w, d_ca(s))) * s.s_ca) .+
                     (Δt*da/2) .* (1 .- s.ϕ) ./ s.ϕ .* (δ .- s.s_ca) .* ξ_3(s),
            s.s_ca - (Δt/(2*Δx)) .* (advection_matrix(w, sigma_peclet(w, d_ca(s))) * s.s_ca) .+
                     (Δt*da/2) .* (1 .- s.ϕ) ./ s.ϕ .* (δ .- s.s_ca) .* ξ_3(s))
        s_co3_half = ifelse.(
            s.ϕ .<= p.ϵ,
            s.s_co3 - (Δt/(2*Δx)) .* (advection_matrix(w, sigma_peclet(w, d_co3(s))) * s.s_co3) .+
                        (Δt*da/2) .* (1 .- s.ϕ) ./ s.ϕ .* (δ .- s.s_co3) .* ξ_3(s),
            s.s_co3 - (Δt/(2*Δx)) .* (advection_matrix(w, sigma_peclet(w, d_co3(s))) * s.s_co3) .+
                        (Δt*da/2) .* (1 .- s.ϕ) ./ s.ϕ .* (δ .- s.s_ca) .* ξ_3(s))
        State(c_a_half, c_c_half, s_ca_half, s_co3_half, ϕ_half)
    end

    half_step
    # pde_c_a(s::State) = -velocity_u(s) * partial(s.c_a) - da .* ξ_1(s)
    # pde_c_c(s::State) = -velocity_u(s) * partial(s.c_c) - da .* ξ_2(s)
    # pde_s_ca(s::State) = 
    #     -velocity_w(s) * partial(s.s_ca) + 
    #     1 ./ s.ϕ * (partial(s.ϕ .* d_ca(s)) * partial(s.s_ca) + s.ϕ .* d_ca(s) * partial_2(s.s_ca)) +
    #     da .* (1 .- s.ϕ) ./ s.ϕ .* (δ .- s.s_ca) .* ξ_3(s)
    # pde_s_co3(s::State) = 
    #     -velocity_w(s) * partial(s.s_co3) + 
    #     1 ./ s.ϕ * (partial(s.ϕ .* d_co3(s)) * partial(s.s_co3) + s.ϕ .* d_co3(s) * partial_2(s.s_co3)) +
    #     da .* (1 .- s.ϕ) ./ s.ϕ .* (δ .- s.s_co3) .* ξ_3(s)
    # pde_ϕ(s::State) =
    #     -partial(velocity_w(s) .* s.ϕ) + d_ϕ(s) * partial_2(s.ϕ) +
    #     da .* (1 - s.ϕ) .* ξ_3(s)
end

function initial_state(p::Param)
    ca_0 = p.s0_ca / sqrt(p.k_c)
    co3_0 = p.s0_co3 / sqrt(p.k_c)
    State(fill(p.c0_a, 1024), fill(p.c0_c, 1024), fill(ca_0, 1024),
          fill(co3_0, 1024), [p.ϕ0; fill(p.ϕ_in, 1023)])
end

function main()
    param = SCENARIO_A(0.5, 0.6)
    half_step = propagators(param, 1e-6, 1.0/1023.0)
    s = initial_state(param)
    println(s)
    s2 = half_step(s)
    println(s2)
end

main()
# ~\~ end
