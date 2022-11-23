---
title: astro-turing-julia
subtitle: some notes
author: Johan Hidding
---

goal: Julia reproduction of L'Heureux paper.

The paper introduces the model, which in the end boils down to four dimensionless equations to be solved:

$$\begin{split}
\partial_t C_A = -U\partial_x C_A - Da\big(&(1 - C_A)C_A (\Omega_{DA} - \nu_1 \Omega_{PA}) \\
        & + \lambda C_A C_C (\Omega_{PC} - \nu_2 \Omega_{DC})\big)
\end{split}$$

$$\begin{split}
\partial_t C_C = -U\partial_x C_C + Da\big(&\lambda(1-C_C)C_C(\Omega_{PC} - \nu_2 \Omega_{DC}) \\ &+ C_A C_C (\Omega_{DA} - \nu_1 \Omega_{PA})\big)
\end{split}$$

$$\begin{split}
\partial_t \hat{c_k} = & -W \partial_x \hat{c_k} + {1 \over \phi} \partial_x\left(\phi d_k \partial_x \hat{c_k}\right) + Da {(1 - \phi) \over \phi} (\delta - \hat{c_k}) \\ & \big(C_A (\Omega_{DA} - \nu_1 \Omega_{PA}) - \lambda C_C (\Omega_{PC} - \nu_2 \Omega_{DC})\big)
\end{split}$$

$$\begin{split}
\partial_t \phi = - \partial_x(W\phi) + d_{\phi} \partial^2_x \phi + Da(1 - \phi) \big(&C_A(\Omega_{DA} - \nu_1\Omega_{PA}) - \\ & \lambda C_C(\Omega_{PC} - \nu_2\Omega_{DC})\big)
\end{split}$$

Where it should be noted that the third equation solves for two quantities, $c_k$ being a rescaled concentration of dissolved species, Ca and CO3. Both these quantities appear in the definitions of $\Omega_i$ (which are therefore functions of $x$). The problem with these equations is that $D_{\phi}$ is very small compared to $D_k$. This causes numerical instabilities when trying to solve these equations using of the shelf PDE methods.

From section 2.8 of the paper we can learn that:

- We should use a Fiadeiro-Veronis scheme to solve the advective diffusion equations, involving $\hat{c_k}$ and $\phi$. The non-linear terms are estimated using the advanced projection method, evaluating terms at $t + \Delta t/2$ by Taylor expansion.
- The advective equations, involving $C_A$ and $C_C$ are solved using an explicit upstream scheme. The gist of that method is that, since things are always moving in the same direction, we may use a forward Euler or RK4 time integration, and a cheap first order finite difference in $x$: the velocity of thy neighbour will soon be yours.
A time step of $\Delta t \sim 5 \times 10^{-6}$ was used.

To see how the Fiadeiro-Veronis scheme works, we follow the book "Diagenetic Models and their Implementation" by Boudreau.

## Crash course in numerical PDE
It is illustrative to break down the very complicated equations we have into smaller parts and see what they do. Let's start with an initial value problem for linear decay:

$$\partial_t C = -k C$$

We know this equation, saying that the rate of change is proportional to the current value of $C$, should result in an exponential solution:

$$C(t) = A \exp(-kt)$$

Setting the initial condition $C(0) = C_0$, we get $A = C_0$. Can we obtain the same through numerical integration?

Starting at $t=0$, we may use a Taylor series expansion to estimate $C$ at $t=\Delta t$:

$$C_1 = C_0 + \partial_t C |_{t=0} \Delta t + O(\Delta t^2).$$

Substituting the ODE into the expansion:

$$C_1 = C_0 - k C_0 \Delta t + O(\Delta t^2).$$

We can generate a series $C_n$:

$$C_{n+1} = C_n (1 - k\Delta t).$$

Now, if we choose $\Delta t$ too big, $\Delta t > 1/k$, then suddenly our generated series will start to oscillate. This is definitely not what we want.

We may go higher order and backward $C(t - \Delta t)$ approximates to:

$$C_i = C_{i+1} - \partial_t C|_{x_{i+1}} \Delta t + \partial_t^2 C|_{x_{i+1}} \Delta t^2 + \dots$$

Leads to,

$$C_{i+1} = C_i + {\Delta t \over 2} \left(\partial_t C|_{x_i} + \partial_t C|_{x_{i+1}}\right),$$

a finite difference scheme very similar to the trapezium rule for function integration. This method is second order implicit.

## Higher order ODE
Suppose we have an equation of the form:

$$D \partial_t^2 C - v \partial_t C + R(C, 0) = 0$$

We need a way to deal with second order derivatives: the central difference formula:

$$C_{i-1} = C_i - \partial C|_i \Delta t + {1\over 2} \partial^2 C|_x \Delta t^2- {1\over 6} \partial^3 C|_i \Delta t^3 + O(\Delta t^4)$$

$$C_{i+1} = C_i + \partial C|_i \Delta t + {1\over 2} \partial^2 C|_x \Delta t^2 + {1\over 6} \partial^3 C|_i \Delta t^3 + O(\Delta t^4)$$

$$C_{i+1} + C_{i-1} = 2C_i + \partial^2 C|_i \Delta t^2 + O(\Delta x^4)$$

$$\partial^2 C|_i = {{C_{i+1} - 2C_i + C_{i-1}} \over \Delta t^2}$$

This leads to a discretisation of our equation as follows,

$$D {{C_{i+1} - 2C_i + C_{i-1}} \over \Delta t^2} - v {{C_{i+1} - C_{i-1}} \over 2 \Delta x} + R(C_i, x_i) = 0$$

Boudreau explains that this is only stable if $v\Delta x / (2D) < 1$. We can make it more stable by using a backwards difference formula, but then we have a very nice second order approximation for the diffusion part and only a first order one for the advection part. Fiadeiro and Veronis (1977) propose the following scheme:

$$D\Big({{C_{i+1} - 2C_{i} + C_{i-1}} \over {\Delta x^2}}\Big) - v \Big({{(1 - \sigma)C_{i+1} + 2\sigma C_i - (1+\sigma)C_{i-1}} \over {2\Delta x}}\Big) = 0,$${#eq:fiadeiro-veronis-scheme}

where, 

$$\sigma = \coth\Big({{v\Delta x} \over {2D}}\Big) - {{2D} \over {v\Delta x}}.$${#eq:fiadeiro-veronis-sigma}

The idea is that $\sigma$ interpolates smoothly between an upstream and central difference method. This behaviour is especially important when, depending on depth (in our 1-d model), the coefficients change: we'll have regions that are dominated by diffusion and others that are dominated by advection.

## Implementation
We take appart the set of equations, and put the complicated non-linear terms (at the ends in brackets) and call them $\xi_1 \dots \xi_3$.

$$\partial_t C_A = -U \partial_x C_A - Da \xi_1,$$
$$\partial_t C_C = -U \partial_x C_C + Da \xi_2,$$
$$\partial_t \hat{c}_k = -W \partial_x \hat{c}_k + {1 \over \phi} \partial_x (\phi d_k \partial_x \hat{c}_k) + Da {{1 - \phi} \over \phi} (\delta - \hat{c}_k) \xi_3,$$
$$\partial_t \phi = -\partial_x(W\phi) + d_{\phi}\partial_x^2\phi + Da(1 - \phi) \xi_3.$$

The first two equations are integrated using an explicit upwind scheme. The non-linear terms in the advection-diffusion equations can be linearised by first-order Taylor approximation, or advanced projection method.

### Parameters
We'll put all the parameters from Table 1 of the model into a big structure `Param` and develop functions on top of that.

::: {.table .table-to-code}
| Parameter             | Scenario A   | Scenario B   | Unit    | Field    |
| --------------------- | ------------ | ------------ | ------- | -------- |
| $\mu_A$               | $100.09$     | $100.09$     | g/mol   | `μ_a`    |
| $\mu_W$               | $18.45$      | $18.45$      | g/mol   | `μ_w`    |
| $\rho_A$              | $2.950$      | $2.950$      | g/cm³   | `ρ_a`    |
| $\rho_C$              | $2.710$      | $2.710$      | g/cm³   | `ρ_c`    |
| $\rho_T$              | $2.8$        | $2.8$        | g/cm³   | `ρ_t`    |
| $\rho_W$              | $1.023$      | $1.023$      | g/cm³   | `ρ_w`    |
| $D_{\rm Ca}^0$        | $131.9$      | $131.9$      | cm²/a   | `d0_ca`  |
| $D_{\rm CO3}^0$       | $272.6$      | $272.6$      | cm²/a   | `d0_co3` |
| $\beta$               | $0.1$        | $0.01$       | cm/a    | `β`      |
| $b$                   | $5$          | $10$         | 1/(kPa) | `b`      |
| $K_A$                 | $10^{-6.19}$ | $10^{-6.19}$ | M²      | `k_a`    |
| $K_C$                 | $10^{-6.37}$ | $10^{-6.37}$ | M²      | `k_c`    |
| $k_1 = k_3$           | $1.0$        | $0.01$       | 1/a     | `k_1`    |
| $k_3 = k_4$           | $0.1$        | $0.001$      | 1/a     | `k_3`    |
| $m = m'$              | $2.48$       | $2.48$       | -       | `m`      |
| $n = n'$              | $2.80$       | $2.80$       | -       | `n`      |
| $x_d$                 | $50$         | $50$         | cm      | `x_d`    |
| $h_d$                 | $100$        | $100$        | cm      | `h_d`    |
| $L$                   | $500$        | $500$        | cm      | `l`      |
| $S$                   | $0.1$        | $0.01$       | cm/a    | `s`      |
| $C_A^0$               | $0.6$        | $0.6$        | -       | `c0_a`   |
| $C_C^0$               | $0.3$        | $0.3$        | -       | `c0_c`   |
| $\hat{c}_{\rm Ca}^0$  | $0.326$      | $0.326$      | mM      | `s0_ca`  |
| $\hat{c}_{\rm CO3}^0$ | $0.326$      | $0.326$      | mM      | `s0_co3` |
:::: {.caption}
**Table 1** Parameters used by L'Heureux 2018.
::::
:::

``` {.julia file=src/model.jl #model}
using Printf: @printf
using LinearAlgebra: Tridiagonal, diagm, I
using HDF5
using ProgressBars

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

SCENARIO_A(ϕ_in :: Float64, ϕ0 :: Float64) = Param(
    100.09, 18.45, 2.950, 2.710, 2.8, 1.023, 131.9, 272.6, 0.1, 5, 10^(-6.19), 10^(-6.37),
    1.0, 1.0, 0.1, 0.1, 2.48, 2.80, 50, 100, 500, 0.1, 0.6, 0.3, 0.326e-3, 0.326e-3, 0.01, 9.81,
    ϕ_in, ϕ0, 0.01)

SCENARIO_B(ϕ_in :: Float64, ϕ0 :: Float64) = Param(
    100.09, 18.45, 2.950, 2.710, 2.8, 1.023, 131.9, 272.6, 0.01, 10, 10^(-6.19), 10^(-6.37),
    0.01, 0.01, 0.001, 0.001, 2.48, 2.80, 50, 100, 500, 0.01, 0.6, 0.3, 0.326e-3, 0.326e-3, 0.01, 9.81,
    ϕ_in, ϕ0, 0.01)

struct State
    c_a :: Vector{Float64}
    c_c :: Vector{Float64}
    s_ca :: Vector{Float64}
    s_co3 :: Vector{Float64}
    ϕ :: Vector{Float64}
end
```

### Numeric difference schemes
In our integration scheme we have essentially three types of discretisation schemes: central difference, Fiadero-Veronis and upwind.

``` {.julia #model}
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

function grad(y::Vector{T}, Δx::T) where T <: Real
    [(y[2] - y[1])/Δx; (y[2:end] - y[1:end-1]) / Δx]
end
```

### Clipping
In the code by L'Heureux there are many cases where values beyond some range, that are considered unphysical, are clipped.

``` {.julia #model}
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
```

### Solver
The largest chunk of code is the actual model. The `propagators` function creates two (closure) functions. One for the half-steps and one for the full-step, that takes the half-step result as an extra argument.

``` {.julia #model}
function propagators(p::Param, Δt::Float64, N::Int)
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

    x = collect(LinRange(0.0, 1.0, N))
    Δx = 1.0 / (N - 1)
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
            s.ϕ + (Δt/(2*Δx))   .* (advection_matrix(1.0, sigma_peclet(w, d_ϕ)) * (s.ϕ .* w)) .+ 
                  (Δt/(2*Δx^2)) .* diffusion_matrix(fill(d_ϕ, N)) * s.ϕ .+
                  (Δt*da/2) .* (1 .- s.ϕ) .* ξ_3(s),
            p.ϵ, 1 - p.ϵ)

        d = d_ca(s)
        vel = w .- grad(d, Δx) - d ./ s.ϕ .* grad(s.ϕ, Δx)
        s_ca_half = s.s_ca + (Δt/(2*Δx)) .* (advection_matrix(vel, sigma_peclet(w, d)) * s.s_ca) .+
                     (Δt/(2*Δx^2)) .* (diffusion_matrix(d ./ s.ϕ) * s.s_ca) .+
                     (Δt*da/2) .* (1 .- s.ϕ) ./ s.ϕ .* (δ .- s.s_ca) .* ξ_3(s)

        d = d_co3(s)
        vel = w .- grad(d, Δx) - d ./ s.ϕ .* grad(s.ϕ, Δx)             
        s_co3_half = s.s_co3 + (Δt/(2*Δx)) .* (advection_matrix(vel, sigma_peclet(w, d)) * s.s_co3) .+
                      (Δt/(2*Δx^2)) .* (diffusion_matrix(d ./ s.ϕ) * s.s_co3) .+
                      (Δt*da/2) .* (1 .- s.ϕ) ./ s.ϕ .* (δ .- s.s_co3) .* ξ_3(s)
        State(c_a_half, c_c_half, s_ca_half, s_co3_half, ϕ_half)
    end

    function full_step(s::State, s_half::State)
        u = velocity_u(s_half)
        w = velocity_w(s_half)

        c_a = s.c_a + (Δt/Δx) .* upwind_dy(s_half.c_a, u) - (da * Δt) .* ξ_1(s_half)
        c_c = s.c_a + (Δt/Δx) .* upwind_dy(s_half.c_c, u) - (da * Δt) .* ξ_2(s_half)
        clip_ara_cal(c_a, c_c)

        _ξ_3 = ξ_3(s_half)
        ϕ_system = I - (Δt/Δx) .* advection_matrix(w, sigma_peclet(w, d_ϕ)) -
                       (Δt/Δx^2) .* diffusion_matrix(fill(d_ϕ, N)) +
                       Δt .* diagm(grad(w, Δx))
        ϕ_react = (Δt*da) .* (1 .- s_half.ϕ) .* _ξ_3
        ϕ_react[1] = 0; ϕ_react[end] = 0
        ϕ = clamp.(ϕ_system\s.ϕ .+ ϕ_react, p.ϵ, 1 - p.ϵ)

        d = d_ca(s_half)
        vel = w .- grad(d, Δx) - x ./ s_half.ϕ .* grad(s_half.ϕ, Δx)
        ca_system = I - (Δt/Δx) .* advection_matrix(vel, sigma_peclet(w, d)) -
                        (Δt/Δx^2) .* diffusion_matrix(d ./ s.ϕ)
        ca_react = (Δt*da) .* (1 .- s_half.ϕ) ./ s_half.ϕ .* (δ .- s_half.s_ca) .* _ξ_3
        ca_react[1] = 0; ca_react[end] = 0
        s_ca = ca_system\s.s_ca .+ ca_react

        d = d_co3(s_half)
        vel = w .- grad(d, Δx) - x ./ s_half.ϕ .* grad(s_half.ϕ, Δx)
        co3_system = I - (Δt/Δx) .* advection_matrix(vel, sigma_peclet(w, d)) -
                         (Δt/Δx^2) .* diffusion_matrix(d ./ s.ϕ)
        co3_react = (Δt*da) .* (1 .- s_half.ϕ) ./ s_half.ϕ .* (δ .- s_half.s_co3) .* _ξ_3
        co3_react[1] = 0; co3_react[end] = 0
        s_co3 = co3_system\s.s_co3 .+ co3_react

        State(c_a, c_c, s_ca, s_co3, ϕ)
    end

    (half_step, full_step)
end

function initial_state(p::Param, N::Int)
    ca_0 = p.s0_ca / sqrt(p.k_c)
    co3_0 = p.s0_co3 / sqrt(p.k_c)
    State(fill(p.c0_a, N), fill(p.c0_c, N), fill(ca_0, N),
          fill(co3_0, N), [p.ϕ0; fill(p.ϕ_in, N-1)])
end

function main()
    fid = h5open("scenario-a.h5", "w")
    param = SCENARIO_A(0.5, 0.6)
    half_step, full_step = propagators(param, 5e-6, 201)
    s = initial_state(param, 201)
    out_aragonite = create_dataset(fid, "aragonite", Float64, (201, 101), chunk=(201,1))
    out_calcite   = create_dataset(fid, "calcite", Float64, (201, 101), chunk=(201,1))
    out_calcium   = create_dataset(fid, "calcium", Float64, (201, 101), chunk=(201,1))
    out_carbonate = create_dataset(fid, "carbonate", Float64, (201, 101), chunk=(201,1))
    out_porosity  = create_dataset(fid, "porosity", Float64, (201, 101), chunk=(201,1))

    function write_snapshot(i::Int, s::State)
        out_aragonite[:,i] = s.c_a
        out_calcite[:,i] = s.c_c
        out_calcium[:,i] = s.s_ca
        out_carbonate[:,i] = s.s_co3
        out_porosity[:,i] = s.ϕ
    end
    write_snapshot(1, s)

    steps = 2000
    save_every = steps ÷ 100
    for i in tqdm(1:steps)
        h = half_step(s)
        s = full_step(s, h)
        if i % save_every == 0
            write_snapshot(i ÷ save_every, s)
        end
    end
    close(fid)
end

main()
```

