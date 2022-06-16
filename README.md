# astro-turing-julia
Julia reproduction of L'Heureux paper.

The paper introduces the model, which in the end boils down to four dimensionless equations to be solved:

$$\partial_t C_A = -U\partial_x C_A - Da\left((1 - C_A)C_A (\Omega_{DA} - \nu_1 \Omega_{PA}) + \lambda C_A C_C (\Omega_{PC} - \nu_2 \Omega_{DC})\right)$$

$$\partial_t C_C = -U\partial_x C_C + Da\left(\lambda(1-C_C)C_C(\Omega_{PC} - \nu_2 \Omega_{DC}) + C_A C_C (\Omega_{DA} - \nu_1 \Omega_{PA})\right)$$

$$\partial_t \hat{c_k} = -W \partial_x \hat{c_k} + {1 \over \phi} \partial_x\left(\phi d_k \partial_x \hat{c_k}\right) + Da {(1 - \phi) \over \phi} (\delta - \hat{c_k}) \left(C_A (\Omega_{DA} - \nu_1 \Omega_{PA}) - \lambda C_C (\Omega_{PC} - \nu_2 \Omega_{DC})\right)$$

$$\partial_t \phi = - \partial_x(W\phi) + d_{\phi} \partial^2_x \phi + Da(1 - \phi) \left(C_A(\Omega_{DA} - \nu_1\Omega_{PA}) - \lambda C_C(\Omega_{PC} - \nu_2\Omega_{DC})\right)$$

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

$$D \partial_t^2 C - v \partial_t C + R(C, 0) = 0$$


