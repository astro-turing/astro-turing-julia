---
title: Advection-diffusion
subtitle: using Fiadeiro-Veronis scheme
---

$$\renewcommand{\vec}[1]{{\bf #1}}$$

# Fiadeiro-Veronis scheme
We start with the ODE,

$$D\partial_x^2 C - v \partial_x C = 0.$${#eq:fv-sample-ode}

Pages 313-314 of Boudreau explains why the previous approaches lead to instabilities here. We need to interpolate between two methods to be stable under a wider set of regimes.

Fiadeiro and Veronis (1977) propose the following scheme:

$${D \over {\Delta x^2}}\Big(C_{i+1} - 2C_{i} + C_{i-1}\Big) - {v \over {2\Delta x}}\Big((1-\sigma)C_{i+1} + 2\sigma C_i - (1+\sigma)C_{i-1}\Big) = 0,$${#eq:fv-scheme}

where,

$$\sigma = \coth\Big({{v\Delta x} \over {2D}}\Big) - {{2D} \over {v \Delta x}}.$${#eq:fv-scheme-sigma}

This means that when we are diffusion dominated ($\sigma = 0$), we have central differencing. When advection becomes dominant ($\sigma = 1$), then we have backwards differencing.

Now, let's have a time-dependent equation,

$$\partial_t x = D\partial_x^2 C - v \partial_x C.$${#eq:fv-sample-pde}

The Fiadeiro-Veronis scheme gives us a tridiagonal matrix when we solve using the implicit Euler method. Similar to when we solved the diffusion equation,

$$C_{j+1} = C_{j} + \Delta t \partial_t C|_{j+1},$${#eq:backward-euler}

and we can write the spatial discretisation as a vector $\vec{C}_i$,

$$\vec{C}_{j+1} = \vec{C}_{j}
  + {{D \Delta t} \over {\Delta x^2}} A_{\rm diff} \vec{C}_{j+1}
  - {{v \Delta t} \over {2 \Delta x}} A_{\rm adv} \vec{C}_{j+1},$${#eq:fv-bwe}

where $A_{\rm diff}$ is the by now familiar $[1, -2, 1]$ tridiagonal matrix, and $A_{\rm adv}$ has a similar structure, as

$$A_{\rm adv} = 
       \begin{pmatrix} 2\sigma &  1-\sigma &  0 &  0 & \dots\\
                       -(1+\sigma) & 2\sigma &  1-\sigma &  0 & \dots\\
                       0 &  -(1+\sigma) & 2\sigma &  1-\sigma & \dots\\
                       0 &  0 &  -(1+\sigma) & 2\sigma & \dots\\
                       \vdots & \vdots & \vdots & \vdots & \ddots\\
\end{pmatrix}$${#eq:adv-matrix}

## Implementation
Let's see how this behaves under variable $D$ and $v$ for our original problem. We define initial conditions with several spikes.

``` {.julia #adv-diff-constants}
const GRID_SIZE = 10000    # grid resolution
const TIME_STEPS = 10000   # number of timesteps
const SIGMA0 = 0.02        # width of peaks
const ADVECTION = 0.26667  # maximum of advection coeff.
const DIFFUSION = 0.01     # maximum of diffusion coeff.
const KAPPA = 0.1          # width of our sigmoid function
```

``` {.julia #adv-diff}
times = LinRange(0.0, 1.0, TIME_STEPS+1)
positions = LinRange(-1.0, 1.0, GRID_SIZE+1)
delta_x = 2.0 / GRID_SIZE
delta_t = 1.0 / TIME_STEPS

y0 = zeros(Float64, GRID_SIZE+1)
for x in -0.8:0.4:0.8
    y0 .+= exp.((positions .- x).^2 ./ (-2*SIGMA0^2))
end
```

We let $D$ and $v$ vary along the interval according to a sigmoid function,

``` {.julia #sigmoid}
function sigmoid(x, kappa=KAPPA)
        1 / (1 + exp(-x / kappa))
end
```

``` {.julia #adv-diff}
diff_coef = DIFFUSION * sigmoid.(positions)
adv_coef = ADVECTION * sigmoid.(-positions)
```

``` {.julia .hide #adv-diff}
@printf "# pos       diffusion           advection\n"
for (x, d, a) in zip(positions, diff_coef, adv_coef)
        @printf "%f %f %f\n" x d a
end
@printf "\n\n"
```

We may rewrite the backwards Euler method to $B C_{j+1} = C_{j}$, where

$$B = I - c_{\rm diff} A_{\rm diff} + c_{\rm adv} A_{\rm adv}.$${#eq:adv-diff-system}

``` {.julia #adv-diff}
c_d = diff_coef .* (delta_t / delta_x^2)
c_a = adv_coef .* (delta_t / (2 * delta_x))
c_f = (adv_coef ./ diff_coef) .* (delta_x / 2)
sigma = coth.(c_f) - 1 ./ c_f
B = Tridiagonal((c_a.*(-sigma .- 1) .- c_d)[2:end],
                1 .+ 2 .* c_d + 2 .* c_a .* sigma,
                (c_a.*(-sigma .+ 1) .- c_d)[1:end-1])
```

We obtain the solution by iteratively solving for $C_{j+1}$.

``` {.julia #adv-diff}
yn = FunctionIterator(y -> B \ y, y0)
result = Iterators.map(Snapshot, times, yn)
```

``` {.julia .hide file=src/adv-diff.jl}
using Printf: @printf
using LinearAlgebra: Tridiagonal

<<adv-diff-constants>>

<<function-iterator>>
<<snapshots>>
<<sigmoid>>

function main()
    <<adv-diff>>
    for s in every_n(result, div(TIME_STEPS, 2))
            print_snapshot(s, positions)
    end
end

main()
```

## Result

``` {.gnuplot .hide file=scripts/adv-diff.gnuplot}
set term svg
set multiplot
set lmargin at screen 0.1
set size 1, 0.35
set yrange [-4:0]
set ylabel "log_{10}(c)"
set xlabel "x"
set key box opaque bottom
plot 'data/adv-diff.out' i 0 u 1:(log10($2)/2) w l lc 6 t'diffusion', \
     '' i 0 u 1:(log10($3)) w l lc 7 t'advection'

set origin 0, 0.32
set size 1, 0.7
unset logscale y
set yrange [-0.1:1.05]
set ylabel "y"
unset xlabel
set label "← advection dominated" at -0.95,-0.05 left
set label "diffusion dominated →" at 0.95,-0.05 right
unset xtics
set key box opaque top
plot 'data/adv-diff.out' i 1 t"t=0.0" w l lc 1, \
     '' i 2 t"t=0.5" w l lc 2, '' i 3 t"t=1.0" w l lc 3
unset multiplot
```

``` {.make .figure target=fig/adv-diff.svg}
Advection-diffusion with Fiadeiro-Veronis scheme.
---
data/adv-diff.out: src/adv-diff.jl
>        @mkdir -p $(@D)
>        julia $< > $@

$(target): scripts/adv-diff.gnuplot data/adv-diff.out
>        gnuplot $< > $@
```

