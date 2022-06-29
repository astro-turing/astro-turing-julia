---
title: Solving the diffusion equation
subtitle: using central differences and backward euler
tags: ["diffusion", "backward euler", "implicit", "pde"]
---

# Some equations
$$\renewcommand{\vec}[1]{{\bf #1}}$$

Let's solve the diffusion equation,

$$\partial_t u = a \partial_{xx} u.$${#eq:diffusion}

We may discretise the double space derivative using central differences,

$$\Delta x^2 \partial_{xx} u|_i = u_{i-1} - 2 u_{i} + u_{i+1}.$${#eq:central-difference}

For the time derivative we use a backward (implicit) Euler scheme,

$$u_{j+1} = u_{j} + \Delta t \partial_t u |_{j+1}.$${#eq:backward-euler}

Combining those gives,

$$u_{(i, j+1)} = u_{(i, j)} + {{a \Delta t} \over {\Delta x^2}} \Big(u_{(i-1, j+1)} - 2 u_{(i, j+1)} + u_{(i+1,j+1)}\Big).$${#eq:diffusion-discretised}

Now, writing our spatially discretized $u_{i}$ as a vector $\vec{u}$, the above relation becomes a matrix equation,

$$\vec{u}_{j+1} = \vec{u}_{j} + {{a \Delta t} \over {\Delta x^2}} A \vec{u}_{j+1},$${#eq:diffusion-matrix-form}

where $A$ is the tri-diagonal matrix,

$$A = \begin{pmatrix} -2 &  1 &  0 &  0 & \dots\\
                      1 & -2 &  1 &  0 & \dots\\
                      0 &  1 & -2 &  1 & \dots\\
                      0 &  0 &  1 & -2 & \dots\\
                      \vdots & \vdots & \vdots & \vdots & \ddots\\
\end{pmatrix}.$${#eq:central-difference-matrix}

This entire system can be written as a linear equation $B \vec{u}_{i+1} = \vec{u}_{i}$, where $B = I - cA$, and $c = a\Delta t / \Delta x^2$. Note that we have assumed boundary conditions where $u = 0$ outside the domain.

# Implementation

``` {.julia file=src/diffusion.jl}
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
```

