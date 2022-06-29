# ~\~ language=Julia filename=src/diffusion.jl
# ~\~ begin <<lit/diffusion.md|src/diffusion.jl>>[init]
using LinearAlgebra: SymTridiagonal
using Printf: @printf

const GRID_SIZE = 50
const TIME_STEPS = 100
const SIGMA0 = 0.1
const DIFFUSION = 0.05

struct FunctionIterator
        f :: Function
        arg :: Any
end

Base.iterate(a::FunctionIterator) = (a.arg, a.f(a.arg))
Base.iterate(a::FunctionIterator, x) = (x, a.f(x))
Base.IteratorSize(::FunctionIterator) = Base.IsInfinite()

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

function every_n(iter, n)
        Iterators.map(first, Iterators.partition(iter, n))
end

function main()
        positions = LinRange(-1.0, 1.0, GRID_SIZE+1)
        times = LinRange(0.0, 1.0, TIME_STEPS+1)
        y0 = exp.(positions.^2 ./ (-2*SIGMA0^2))
        delta_x = 2.0 / GRID_SIZE
        delta_t = 1.0 / TIME_STEPS
        c = DIFFUSION * delta_t / delta_x^2 
        B = SymTridiagonal(fill(1 + 2*c, GRID_SIZE+1), fill(-c, GRID_SIZE))

        result = Iterators.map(Snapshot, times, FunctionIterator(y -> B \ y, y0))
        for s in every_n(result, 50)
                print_snapshot(s, positions)
        end
end

main()
# ~\~ end
