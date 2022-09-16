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
We take appart the set of equations, and put the complicated non-linear terms (at the ends in brackets) and call them $\xi_1 \dots \xi_4$. The first two equations are integrated using an explicit upwind scheme.

$$\partial_t C_A = -U \partial_x C_A - Da \xi_1,$$
$$\partial_t C_C = -U \partial_x C_C + Da \xi_2,$$
$$\partial_t \hat{c}_k = -W \partial_x \hat{c}_k + {1 \over \phi} \partial_x (\phi d_k \partial_x \hat{c}_k) + Da {{1 - \phi} \over \phi} (\delta - \hat{c}_k) \xi_3,$$
$$\partial_t \phi = -\partial_x(W\phi) + d_{\phi}\partial_x^2\phi + Da(1 - \phi) \xi_4.$$

We'll put all the parameters from Table 1 of the model into a big structure `Param` and develop functions on top of that.

::: {.table .table-to-code}
| Parameter             | Scenario A   | Scenario B   | Unit    | Field    |
| ---------             | ----------   | ----------   | ----    | -----    |
| $\mu_A$               | $100.09$     | $100.09$     | g/mol   | `mu_a`   |
| $\mu_W$               | $18.45$      | $18.45$      | g/mol   | `mu_w`   |
| $\rho_A$              | $2.950$      | $2.950$      | g/cm³   | `rho_a`  |
| $\rho_C$              | $2.710$      | $2.710$      | g/cm³   | `rho_c`  |
| $\rho_T$              | $2.8$        | $2.8$        | g/cm³   | `rho_t`  |
| $D_{\rm Ca}^0$        | $131.9$      | $131.9$      | cm²/a   | `d0_ca`  |
| $D_{\rm CO3}^0$       | $272.6$      | $272.6$      | cm²/a   | `d0_co3` |
| $\beta$               | $0.1$        | $0.01$       | cm/a    | `beta`   |
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

``` {.julia file=src/model.jl}
struct Param
    mu_a :: Float64
    mu_w :: Float64
    rho_a :: Float64
    rho_c :: Float64
    rho_t :: Float64
    rho_w :: Float64
    d0_ca :: Float64
    d0_co3 :: Float64
    beta :: Float64
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
    phi0 :: Float64
    phi_in :: Float64            # initial porosity
    phi_inf :: Float64           # set to 0.01, see p. 7 bottom left
    g :: Float64                 # gravity, set to 9.8 m/s^2
end

struct Aux
    x :: Vector{Float64}
    theta :: Vector{Float64}
end

struct State
    c_a :: Vector{Float64}
    c_c :: Vector{Float64}
    s_ca :: Vector{Float64}
    s_co3 :: Vector{Float64}
    phi :: Vector{Float64}
end

damkoehler_number(p::Param) = p.k_2 * p.d0_ca / p.s^2
lambda(p::Param) = p.k_3 / p.k_2
nu_1(p::Param) = p.k_1 / p.k_2
nu_2(p::Param) = p.k_4 / p.k_3
d_ca(p::Param, s::State) = 1.0 ./ (1 .- 2 .* log.(s.phi))
d_co3(p::Param, s::State) = (p.d0_co3 / p.d0_ca) ./ (1 .- 2 .* log.(s.phi))

# Equations 15: K
conductivity(p::Param, phi::Float64) = beta * phi^3 / (1 - phi)^2 * conductivity_correction(phi)

# Equation 17: F
conductivity_correction(phi::Float64) = 1 - exp(10 - 10/phi)

# Equation 24: H
elastic_response(p::Param, phi_nr::Float64) = -1 / (p.b * (phi_nr - p.phi_inf))

# Equation 25: D_phi, this function is taken as a constant in the diffusion terms
# and is completely ignored in the advection terms.
porosity_diffusion(p::Param, phi::Float64) =
    (p.beta / p.b * p.g * p.rho_w) * phi^3 / (1 - phi) / (p.phi_in - p.phi_inf) *
    conductivity_correction(phi)

# This is strictly speaking a function of phi, but taken to be constant
d_phi(p::Param) = porosity_diffusion(p, p.phi_in) / p.d0_ca
x_star(p::Param) = p.d0_ca / p.s
t_star(p::Param) = p.d0_ca / p.s^2
function theta(p::Param, x::Float64)
    x_prime = x * x_star(p)
    if p.x_d < x_prime && x_prime < (p.x_d + p.h_d)
        1
    else
        0
    end
end

# § 2.5 states that rho_s is taken as a constant, therefore delta will be a constant
rho_s(p::Param) = p.rho_t / (1 - p.c0_a * (1 - p.rho_t / p.rho_a) - p.c0_c * (1 - p.rho_t / p.rho_c))
delta(p::Param) = rho_s(p) / (p.mu_a * sqrt(p.k_c))

# Equation 45: TODO These contain many common subexpressions that can be eliminated
omega_pa(p::Param, s::State) = (s.s_ca .* s.s_co3 .* (p.k_c / p.k_a) .- 1).^(p.m)
omega_da(p::Param, s::State, aux::Aux) = (1 .- s.s_ca .* s.s_co3 .* (p.k_c / p.k_a)).^(p.m) .* aux.theta
omega_pc(p::Param, s::State) = (s.s_ca .* s.s_co3 .- 1).^(p.n)
omega_dc(p::Param, s::State) = (1 .- s.s_ca .* s.s_co3).^(p.n)

# Equations 46 & 47: TODO reduce common subexpressions
velocity_u(p::Param, s::State) =
    1 - conductivity(p, p.phi0) / p.s * (1 - p.phi0) * (rho_s(p) / p.rho_w - 1) .+
    conductivity.(p, s.phi) ./ p.s .* (1 .- s.phi) .* (rho_s(p) / p.rho_w - 1)

velocity_w(p::Param, s::State) =
    1 - conductivity(p, p.phi0) / p.s * (1 - p.phi0) * (rho_s(p) / p.rho_w - 1) .+
    conductivity.(p, s.phi) ./ p.s .* (1 .- s.phi)^2 ./ s.phi .* (rho_s(p) / p.rho_w - 1)

function partial(y)
    y   # TODO implement some representation of a differential operator
end

xi_1(p::Param, s::State, aux::Aux) =
    (1 .- s.c_a) .* s.c_a .* (omega_da(p, s, a) .- nu_1(p) .* omega_pa(p, s)) +
    lambda(p) .* s.c_a .* s.c_c .* (omega_pc(p, s) - nu_2(p) .* omega_dc(p, s))

xi_2(p::Param, s::State, aux::Aux) = 

pde_c_a(p::Param, s::State, aux::Aux) = -velocity_u(p, s) * partial(s.c_a) - damkoehler_number(p) .* xi_1(p, s, aux)
pde_c_c(p::Param, s::State, aux::Aux) = -velocity_u(p, s) * partial(s.c_c) - damkoehler_number(p) .* xi_2(p, s, aux)

# TODO expand derivative => partial(s.phi .* d_k(p)) * partial(s.s_<k>) + s.phi .* d_k(p) * partial_2(s.s_<k>)
pde_s_ca(p::Param, s::State, aux::Aux) = 
    -velocity_w(p, s) * partial(s.s_ca) + 
    1 ./ s.phi * partial(s.phi .* d_ca(p, s) * partial(s.s_ca)) +
    damkoehler_number(p) .* (1 .- s.phi) ./ s.phi .* (delta(p) .- s.s_ca) .* xi_3(p, s, aux)

pde_s_co3(p::Param, s::State, aux::Aux) = 
    -velocity_w(p, s) * partial(s.s_co3) + 
    1 ./ s.phi * partial(s.phi .* d_co3(p, s) * partial(s.s_co3)) +
    damkoehler_number(p) .* (1 .- s.phi) ./ s.phi .* (delta(p) .- s.s_co3) .* xi_3(p, s, aux)

pde_phi(p::Param, s::State, aux::Aux) =
    -partial(velocity_w(p, s) .* s.phi) + d_phi(p) * partial_2(s.phi) +
    damkoehler_number(p) .* (1 - s.phi) .* xi_4(p, s, aux)

function main()
end

main()
```

