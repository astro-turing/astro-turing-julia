# ~\~ language=Julia filename=src/adv-diff.jl
# ~\~ begin <<lit/fiadeiro-veronis.md|src/adv-diff.jl>>[init]
using Printf: @printf
using LinearAlgebra: Tridiagonal

# ~\~ begin <<lit/fiadeiro-veronis.md|adv-diff-constants>>[init]
const GRID_SIZE = 10000    # grid resolution
const TIME_STEPS = 10000   # number of timesteps
const SIGMA0 = 0.02        # width of peaks
const ADVECTION = 0.26667  # maximum of advection coeff.
const DIFFUSION = 0.01     # maximum of diffusion coeff.
const KAPPA = 0.1          # width of our sigmoid function
# ~\~ end

# ~\~ begin <<lit/diffusion.md|function-iterator>>[init]
struct FunctionIterator
        f :: Function
        arg :: Any
end

Base.iterate(a::FunctionIterator) = (a.arg, a.f(a.arg))
Base.iterate(a::FunctionIterator, x) = (x, a.f(x))
Base.IteratorSize(::FunctionIterator) = Base.IsInfinite()

function every_n(iter, n)
        Iterators.map(first, Iterators.partition(iter, n))
end
# ~\~ end
# ~\~ begin <<lit/diffusion.md|snapshots>>[init]
struct Snapshot
        t :: Float64
        y :: Vector{Float64}
end

function print_snapshot(s::Snapshot, pos)
        @printf "# t = %f\n" s.t
        for (x, y) in zip(pos, s.y)
                @printf "%f %f\n" x y
        end
        @printf "\n\n"
end
# ~\~ end
# ~\~ begin <<lit/fiadeiro-veronis.md|sigmoid>>[init]
function sigmoid(x, kappa=KAPPA)
        1 / (1 + exp(-x / kappa))
end
# ~\~ end

function main()
    # ~\~ begin <<lit/fiadeiro-veronis.md|adv-diff>>[init]
    times = LinRange(0.0, 1.0, TIME_STEPS+1)
    positions = LinRange(-1.0, 1.0, GRID_SIZE+1)
    delta_x = 2.0 / GRID_SIZE
    delta_t = 1.0 / TIME_STEPS

    y0 = zeros(Float64, GRID_SIZE+1)
    for x in -0.8:0.4:0.8
        y0 .+= exp.((positions .- x).^2 ./ (-2*SIGMA0^2))
    end
    # ~\~ end
    # ~\~ begin <<lit/fiadeiro-veronis.md|adv-diff>>[1]
    diff_coef = DIFFUSION * sigmoid.(positions)
    adv_coef = ADVECTION * sigmoid.(-positions)
    # ~\~ end
    # ~\~ begin <<lit/fiadeiro-veronis.md|adv-diff>>[2]
    @printf "# pos       diffusion           advection\n"
    for (x, d, a) in zip(positions, diff_coef, adv_coef)
            @printf "%f %f %f\n" x d a
    end
    @printf "\n\n"
    # ~\~ end
    # ~\~ begin <<lit/fiadeiro-veronis.md|adv-diff>>[3]
    c_d = diff_coef .* (delta_t / delta_x^2)
    c_a = adv_coef .* (delta_t / (2 * delta_x))
    c_f = (adv_coef ./ diff_coef) .* (delta_x / 2)
    sigma = coth.(c_f) - 1 ./ c_f
    B = Tridiagonal((c_a.*(-sigma .- 1) .- c_d)[2:end],
                    1 .+ 2 .* c_d + 2 .* c_a .* sigma,
                    (c_a.*(-sigma .+ 1) .- c_d)[1:end-1])
    # ~\~ end
    # ~\~ begin <<lit/fiadeiro-veronis.md|adv-diff>>[4]
    yn = FunctionIterator(y -> B \ y, y0)
    result = Iterators.map(Snapshot, times, yn)
    # ~\~ end
    for s in every_n(result, div(TIME_STEPS, 2))
            print_snapshot(s, positions)
    end
end

main()
# ~\~ end
