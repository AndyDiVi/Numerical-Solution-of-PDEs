# Numerical Treatment of Partial Differential Equations

Authors: Andrea Di Via, Keylor Bazzotti & Michele Lai  


https://github.com/user-attachments/assets/493bf24c-57d6-4072-a403-0a04aed185b9



## Overview

This project studies the numerical solution of elliptic and parabolic partial differential equations (PDEs) on a two-dimensional domain with mixed boundary conditions.

The work covers:

- The Laplace equation for a stationary problem without internal sources.
- A Poisson-type problem with localized sources and piecewise material coefficients.
- The heat equation for time-dependent diffusion.
- Finite Difference Method (FDM) implementations in MATLAB.
- Finite Element Method (FEM) comparisons using MATLAB's PDE Toolbox.
- Sparse linear systems, boundary treatment, and mesh-resolution experiments.

This README summarizes the accompanying report. Script names and repository layout are not specified in the report, so the execution instructions do not assume particular filenames.

## Domain and Geometry

The domain is symmetric about both coordinate axes. Its geometry is described in the report using:

```math
P=(12,8), \qquad Q=(4,2), \qquad R=(6,2)
```

The notation distinguishes:

- `Omega_1`: the surrounding region.
- `Omega_2`: the internal source/material regions.
- `Gamma_int`: internal boundaries in the Laplace problem.
- `Gamma_ext`: the external boundary, divided into north, south, east, and west sides.

The Laplace problem is solved on the surrounding region. The Poisson and heat problems include the internal regions as part of the computational domain.

Refer to the domain diagrams in the report for the full boundary geometry and interface layout.

## Mathematical Models

### 1. Laplace Equation

The stationary source-free problem is:

```math
\Delta u = 0 \qquad \text{in }\Omega_1
```

The north and south boundaries have prescribed values:

```math
u=g_0(x), \qquad g_0(x)=10(12-|x|)^2
```

The east boundary has:

```math
u=4
```

The west boundary has a homogeneous Neumann condition:

```math
\nabla u\cdot n=0
```

The internal boundaries have:

```math
\nabla u\cdot n=100
```

Here, the normal direction follows the boundary convention of the problem. At connecting vertices, the report gives priority to the north/south Dirichlet data.

### 2. Poisson-Type Problem

The report denotes the stationary source problem by:

```math
-K\Delta u=f
```

The source is localized in the internal regions:

```math
f(x,y)=\begin{cases}0 & (x,y)\in\Omega_1\\1000 & (x,y)\in\Omega_2\end{cases}
```

The coefficient is:

```math
K(x,y)=\begin{cases}1 & (x,y)\in\Omega_1\\K_0 & (x,y)\in\Omega_2,\ x>0\\K_0^{-1} & (x,y)\in\Omega_2,\ x\leq0\end{cases}
```

Two cases are compared:

- `K0 = 1`: constant coefficient.
- `K0 = 110`: discontinuous coefficient with high material contrast.

The external boundary conditions are the same as in the Laplace problem: Dirichlet data on north, south, and east, and homogeneous Neumann data on west. The former internal boundaries become interfaces.

For the discontinuous-coefficient case, the report introduces coefficient-dependent and interface stencils. Reproducing the results requires following those stencils, not simply multiplying every row of a constant-coefficient Laplacian by `K`. The distinction between the displayed non-divergence form and a conservative variable-conductivity model should be checked when extending the implementation.

### 3. Heat Equation

The transient problem is:

```math
u_t-\frac{1}{100}\Delta u=0
```

The initial condition is the finite-difference solution of the constant-coefficient Poisson problem:

```math
u(x,y,0)=u_0(x,y)
```

The external boundary conditions remain fixed for positive times.

The main experiments use the time interval `[0, 10]`. The report also discusses a longer simulation up to `T = 600` to make the slow diffusion more visible.

## Numerical Methods

### Finite Differences

A uniform Cartesian grid is constructed with spatial step `h`. The reported experiments use:

```math
h\in\{1,0.5,0.2\}
```

The interior Laplacian is approximated using the five-point stencil:

```math
\Delta_hu_{i,j}=\frac{u_{i+1,j}+u_{i-1,j}+u_{i,j+1}+u_{i,j-1}-4u_{i,j}}{h^2}
```

The implementation uses two auxiliary matrices:

- A geometry matrix classifying grid points as interior, boundary, interface, outside-domain, or ghost points.
- An index matrix mapping unknown grid values to equations in the linear system.

Known Dirichlet values contribute to the right-hand side. Neumann conditions are treated using centered differences and ghost points.

The stationary discretization gives:

```math
A_hU=b
```

The matrix is assembled in sparse format and solved using MATLAB's backslash operator:

```matlab
U = Ah \ b;
```

This solves the linear system without explicitly forming the inverse matrix.

### Finite Elements

The elliptic problems are also solved using triangular meshes through MATLAB's `pdetool` interface.

The FEM workflow includes:

1. Defining the polygonal geometry.
2. Assigning boundary conditions.
3. Setting PDE coefficients and sources.
4. Generating and refining the mesh.
5. Solving the discrete problem.
6. Evaluating the FEM solution at selected comparison points.

The report uses piecewise-linear finite element approximations and describes interpolation through `pdeInterpolant` and `evaluate`.

### Crank–Nicolson Time Integration

The heat equation is primarily solved with Crank–Nicolson, using:

```math
k=\frac{h}{2}
```

The diffusion ratio is:

```math
\lambda=\frac{c^2k}{h^2}, \qquad c^2=\frac{1}{100}
```

The discrete system at each time step is:

```math
A_hU^{n+1}=B_hU^n+b
```

For interior points, the left-hand stencil has center coefficient `1 + 2*lambda` and neighbor coefficients `-lambda/2`. The right-hand stencil has center coefficient `1 - 2*lambda` and neighbor coefficients `lambda/2`.

The formal interior discretization error is:

```math
O(k^2)+O(h^2)
```

Boundary treatment, interface regularity, and geometry must also be considered when assessing overall convergence.

Explicit Euler and implicit Euler are implemented as additional comparisons for the transient problem. The report states that their outputs showed no significant differences in the tested configurations.

## Reported Results

### Laplace Problem

The reported solution has:

- Maximum value `1440` at `(0, ±8)`.
- Minimum value `0` at the north/south corner points `(±12, ±8)`.

The report also analyzes the negative gradient:

```math
E=-\nabla u
```

#### Linear System Size

| Spatial step | Matrix dimensions |
|---|---:|
| 1 | 375 × 375 |
| 0.5 | 1509 × 1509 |
| 0.2 | 9321 × 9321 |

The reported FDM–FEM discrepancies do not decrease consistently under mesh refinement for this case. These results should be treated as a numerical comparison requiring further verification, rather than evidence of established convergence.

### Poisson Problem

The maximum occurs at `(-5, 0)` in the reported experiments.

| Spatial step | Maximum for K0 = 110 | Maximum for K0 = 1 |
|---|---:|---:|
| 1 | Approximately 49400 | Approximately 2830 |
| 0.5 | Approximately 52300 | Approximately 3450 |
| 0.2 | Approximately 53400 | Approximately 3740 |

The low-coefficient internal region exhibits a much larger temperature peak in the discontinuous case.

The reported pointwise FDM–FEM differences generally decrease as the grid is refined for the Poisson experiments. These differences compare two numerical approximations; neither is an exact reference solution.

### Heat Equation

The maximum temperature decreases gradually over the main simulation interval:

| Spatial step | Maximum at t = 0 | Maximum at t = 10 |
|---|---:|---:|
| 1 | Approximately 2830 | Approximately 2790 |
| 0.5 | Approximately 3450 | Approximately 3420 |
| 0.2 | Approximately 3740 | Approximately 3640 |

The relatively small diffusivity makes evolution over `[0, 10]` slow. Longer simulations in the report illustrate the subsequent redistribution of temperature.

## Visualizations

The report includes:

- Domain geometry and boundary classification.
- Cartesian grids and FEM triangulations.
- Geometry matrices and sparse matrix structure.
- Two-dimensional solution maps.
- Three-dimensional solution surfaces.
- Negative-gradient vector fields.
- Pointwise FDM–FEM comparison tables.
- Time-dependent heat diffusion results.

## Requirements

- MATLAB for the finite-difference computations and visualizations.
- Partial Differential Equation Toolbox for the FEM workflow using `pdetool`.
- The accompanying project scripts and geometry definitions.

The report documents the methods but does not provide a complete repository file inventory. Verify the available scripts before running the experiments.

## How to Run

1. Open the project folder in MATLAB.
2. Identify the scripts corresponding to Laplace, Poisson, and heat-equation experiments.
3. Set the domain geometry and spatial step.
4. Run the stationary finite-difference experiments.
5. For Poisson, compare `K0 = 1` and `K0 = 110`.
6. Use the constant-coefficient Poisson solution as the initial condition for the heat equation.
7. Set the final time and use `k = h/2` for the documented Crank–Nicolson setup.
8. Reproduce the elliptic FEM comparisons in PDE Toolbox using matching geometries and boundary conventions.
9. Compare solution values, extrema, and plots across resolutions.

Actual script names depend on the files supplied with the repository.

## Interpretation and Limitations

- The experiments use numerical approximations, not closed-form reference solutions.
- FDM–FEM agreement alone does not prove correctness or convergence.
- The Laplace comparisons exhibit unresolved mesh-dependent discrepancies.
- Discontinuous coefficients require careful interface treatment and a consistent PDE formulation.
- FDM grid spacing and FEM element size are different mesh descriptors and do not imply identical degrees of freedom.
- Reported timings depend on hardware, MATLAB version, and implementation details.
- Overflowing determinant values are not a reliable numerical test of nonsingularity.

## Conclusions

The project combines sparse finite-difference discretization, finite-element comparisons, and implicit time integration to study stationary and transient PDEs with mixed boundary conditions.

The Poisson experiments illustrate the effect of localized sources and strong material contrast. The heat-equation experiments show slow diffusion from a stationary source-generated initial state after the sources are removed.

The results also highlight the importance of boundary and interface implementation, mesh-resolution studies, and cautious interpretation of comparisons between numerical methods.

---

Numerical Treatment of Differential Equations — Group 01  
Andrea Di Via, Keylor Bazzotti & Michele Lai — Università di Genova, 2023/2024
