# ~\~ language=Julia filename=src/diffusion.jl
# ~\~ begin <<lit/diffusion.md|src/diffusion.jl>>[init]
using LinearAlgebra: SymTridiagonal
using Printf: @printf

# ~\~ begin <<lit/diffusion.md|diffusion-constants>>[init]
const GRID_SIZE = 150
const TIME_STEPS = 100
const SIGMA0 = 0.1
const DIFFUSION = 0.05
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

function main()
        # ~\~ begin <<lit/diffusion.md|diffusion-linear-system>>[init]
        positions = LinRange(-1.0, 1.0, GRID_SIZE+1)
        times = LinRange(0.0, 1.0, TIME_STEPS+1)
        delta_x = 2.0 / GRID_SIZE
        delta_t = 1.0 / TIME_STEPS

        y0 = exp.(positions.^2 ./ (-2*SIGMA0^2))
        # y0 = positions .|> x -> abs(x) > SIGMA0 ? 0 : 1   # alternative: step-function
        c = DIFFUSION * delta_t / delta_x^2 
        B = SymTridiagonal(fill(1 + 2*c, GRID_SIZE+1), fill(-c, GRID_SIZE))
        yn = FunctionIterator(y -> B \ y, y0)
        # ~\~ end
        result = Iterators.map(Snapshot, times, yn)

        for s in every_n(result, 50)
                print_snapshot(s, positions)
        end
end

main()
# ~\~ end
