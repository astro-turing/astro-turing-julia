---
title: "Julia implementation of some numerical PDE solvers"
subtitle: "diagenetic modelling in Julia"
include-after: "[Netherlands eScience Center](https://esciencecenter.nl) and [Utrecht University, faculty of geosciences](https://www.uu.nl/organisatie/faculteit-geowetenschappen)"
---

[![Entangled badge](https://img.shields.io/badge/entangled-Use%20the%20source!-%2300aeff)](https://entangled.github.io/)

# Diagenetic modelling in Julia
This repository contains example codes to solve numerical PDEs in Julia. The goal is to reproduce the results from Ivan l'Heureux 2018. We'll work from basics of solving PDEs to running the full model.

- [Diffusion equation](./diffusion.html)
- [Advection equation](./upwind-scheme.html)
- [Fiadeiro-Veronis scheme](./fiadeiro-veronis.html)
- [l'Heureux model](./lheureux-model.html)

# Architecture
In this implementation, readability is a primary concern. The target model has five long equations and more than 30 free parameters. It would be a shame if we have to mix the intricacies of numerical methods with the actual implementation of the model. When solving partial differential equations, this means that we may want to write the actual PDEs in a natural form, and then automatically generate the space descretisation for each equation. There is a set of Julia package called `ModelingToolkit.jl` and more specifically `MethodOfLines.jl` that implements a scheme like this, so we could end up using that. For now, using those would fall outside of the scope of the project due to time constraints (Learning how to use these and adapting the model could take a significant effort).

The `ModelingToolkit` has a nice [tutorial on "Large Stiff Equations"](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/advanced_ode_example/) that could be of use here.
