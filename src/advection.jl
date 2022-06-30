# ~\~ language=Julia filename=src/advection.jl
# ~\~ begin <<lit/upwind-scheme.md|src/advection.jl>>[init]
using Printf: @printf

# ~\~ begin <<lit/upwind-scheme.md|advection-constants>>[init]
const GRID_SIZE = 50
const TIME_STEPS = 100
const SIGMA0 = 0.1
const ADVECTION = 1.0
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
# ~\~ begin <<lit/upwind-scheme.md|upwind-delta-y>>[init]
function delta_y(y::Array{T}, a) where T
        d = Array{T}(undef, length(y))
        for (i, a) in enumerate(a)
                if i == 1 || i == length(y)
                        d[i] = 0
                else
                        d[i] = a < 0 ? y[i+1] - y[i] : y[i] - y[i-1]
                end
        end
        d
end
# ~\~ end

function main()
        # ~\~ begin <<lit/upwind-scheme.md|advection-upwind-scheme>>[init]
        positions = LinRange(-1.0, 1.0, GRID_SIZE+1)
        times = LinRange(0.0, 1.0, TIME_STEPS+1)
        delta_x = 2.0 / GRID_SIZE
        delta_t = 1.0 / TIME_STEPS

        y0 = exp.((positions .- 0.2).^2 ./ (-2*SIGMA0^2))
        # y0 = positions .|> x -> abs(x - 0.2) > SIGMA0 ? 0.0 : 1.0   # alternative: step-function
        c = ADVECTION * sin.(positions .* pi) * delta_t / delta_x
        f = y -> y .- c .* delta_y(y, c)
        yn = FunctionIterator(f, y0)
        # ~\~ end
        result = Iterators.map(Snapshot, times, yn)
        for s in every_n(result, 50)
                print_snapshot(s, positions)
        end
end

main()
# ~\~ end
