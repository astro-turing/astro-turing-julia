# ~\~ language=Julia filename=src/model.jl
# ~\~ begin <<lit/lheureux-model.md|src/model.jl>>[init]
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
# ~\~ end
