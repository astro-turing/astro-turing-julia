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

``` {.julia file=src/advection-diffusion.jl}

```
