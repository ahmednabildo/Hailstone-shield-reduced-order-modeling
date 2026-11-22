# POD-Based Reduced-Order Modeling for Hailstone Shield Design

**Finite element simulation, Proper Orthogonal Decomposition (POD), and reduced-order design optimization of a two-layer protective shield subjected to a simplified hailstone load.**

This project investigates how **projection-based model order reduction** can replace repeated, relatively expensive finite element solves with a fast, accurate reduced-order model (ROM), following classical POD and reduced-basis concepts [4, 5]. The ROM is then used to find a minimum top-layer stiffness that satisfies a maximum vertical-deformation constraint over a range of impact angles.

> **Scope:** This is a **static linear-elasticity** model with an applied boundary load representing hailstone impact. It does **not** simulate transient contact, impact dynamics, damage, or fracture.

## Key results

The results below are taken from the executed project notebook [1].

| Metric | Recorded notebook result |
| --- | ---: |
| Full-order simulations (snapshots) | 100 |
| Training / test snapshots | 50 / 50 |
| Retained POD modes | 10 |
| Energy captured by 10 modes | 99.99999999999243% |
| Full-order solve time | 148.98 ms |
| Reduced-order solve time | 0.89 ms |
| Per-query speedup | **~167×** |
| Mean relative error on the test snapshots | **4.74295 × 10⁻⁷** |

*These figures are outputs saved in the notebook, not independently reproduced benchmarks. Timing covers one FOM or ROM evaluation and excludes offline snapshot generation, basis construction, and operator projection.*

## Problem setup

The shield is modeled as a **2D, two-layer square domain** (originally 6 cm × 6 cm, rescaled to the unit square). The bottom edge is clamped; the top edge receives a localized, angle-dependent load. Gravity acts throughout the domain. The problem setup and implementation are documented in the assignment statement and notebook [1, 2].

- **Geometry:** unit square, split into upper and lower material layers at `y = 0.5`
- **Finite elements:** 40 × 40 mesh; vector-valued, first-order continuous Galerkin displacement space
- **Lower layer:** fixed reference stiffness
- **Upper layer:** variable stiffness factor `r`
- **Impact region:** centered at `x = 0.5` on the top edge, with normalized half-width `0.1`
- **Boundary condition:** zero displacement on the bottom edge
- **Loading:** constant downward body force plus localized horizontal and vertical impact-load components

The parameter vector is

$$
\boldsymbol{\mu}=(r,m,\theta),
$$

where `r` scales upper-layer stiffness, `m` is the impact-load magnitude, and `θ` is the impact angle in radians.

| Parameter | Sampling range | Meaning |
| --- | --- | --- |
| `r` | 0.5 to 1.0 | Upper-layer stiffness factor |
| `m` | 0.5 to 2.0 | Applied impact-load magnitude |
| `θ` | −0.4 to 0.4 rad | Impact angle |

## Methodology

### 1. Full-order finite element model

The notebook assembles parameter-independent finite element operators using **legacy FEniCS** [3] and exploits their **affine decomposition**, a standard building block of efficient reduced-basis methods [5]:

$$
A(\boldsymbol{\mu}) = rA_1 + A_2
$$

$$
F(\boldsymbol{\mu}) = F_1 + m\sin(\theta)F_2 - m\cos(\theta)F_3
$$

Here, `A₁` and `A₂` represent the upper- and lower-layer contributions to the stiffness matrix, respectively. The load vectors `F₁`, `F₂`, and `F₃` represent gravity and the horizontal/vertical shapes of the localized impact load.

The full-order displacement vector is computed from the constrained system:

$$
A(\boldsymbol{\mu})u(\boldsymbol{\mu})=F(\boldsymbol{\mu}).
$$

The notebook's `Affine()` routine preassembles the operators; `FOMsolver()` applies the bottom-edge Dirichlet condition and solves the resulting sparse linear system.

### 2. Snapshot generation and POD

With `np.random.seed(0)`, the notebook draws **100** parameter combinations from the ranges above and computes the corresponding full-order displacement snapshots. The first **50** snapshots are used to construct the POD basis; the remaining **50** are held out for testing.

Singular value decomposition (SVD) of the training snapshot matrix provides the dominant spatial modes, following the snapshot-based POD methodology [4]:

$$
S_{\mathrm{train}} = U\Sigma W^T, \qquad V=U_{[:,1:N]}, \quad N=10.
$$

The fraction of snapshot energy retained by the basis is

$$
E_N=\frac{\sum_{i=1}^{N}\sigma_i^2}{\sum_i\sigma_i^2}.
$$

The notebook reports **99.99999999999243%** captured energy with ten modes.

### 3. Reduced-order model: offline–online decomposition

The full-order operators are projected onto the POD space **once** during the offline stage. This is the reduced-basis offline–online decomposition discussed in [5]:

$$
A_{i,\mathrm{rb}}=V^T A_iV,
\qquad
F_{j,\mathrm{rb}}=V^TF_j.
$$

For any new parameter vector, the online stage assembles and solves only the small reduced system:

$$
A_{\mathrm{rb}}(\boldsymbol{\mu})c(\boldsymbol{\mu})
=F_{\mathrm{rb}}(\boldsymbol{\mu}),
\qquad
u_{\mathrm{rb}}(\boldsymbol{\mu})=Vc(\boldsymbol{\mu}).
$$

This avoids assembling and solving the full finite element system for every new set of impact conditions.

### 4. ROM validation and speed comparison

The notebook compares full-order and reduced-order displacement fields, benchmarks a single evaluation of each solver, and computes mean relative error across the held-out test set:

$$
\mathrm{MRE}=\frac{1}{N_{\mathrm{test}}}
\sum_{k=1}^{N_{\mathrm{test}}}
\frac{\|u_k-u_{\mathrm{rb},k}\|_2}{\|u_k\|_2}.
$$

The saved result in [1] is **4.742947907880425 × 10⁻⁷**, with measured solve times of **148.98 ms** (FOM) and **0.89 ms** (ROM), approximately a **167×** per-evaluation speedup.

### 5. Stiffness optimization under an impact-angle constraint

The final task uses the ROM to search for the **smallest upper-layer stiffness** satisfying a critical maximum vertical-displacement requirement [1, 2]:

$$
\min_r\;r
\quad\text{subject to}\quad
\max_{(x,y)\in\Omega}|u_y(x,y;r,m_0,\theta)|\leq\Delta_{\mathrm{crit}}.
$$

The notebook uses:

- Fixed impact-load magnitude: `m₀ = 1.5`
- Critical vertical-deformation threshold: `Δcrit = 0.5`
- **100** sampled angles between `−0.4` and `0.4` radians
- Stiffness search interval: `r ∈ [0.5, 1.0]`
- Bisection root-finding with `xtol = 1e-4`

For each angle, the bisection solver searches for the stiffness at which maximum vertical displacement reaches the limit. The final plot shows the required stiffness versus impact angle and a **sampled worst-case design** obtained from the largest required stiffness on the angle grid.

The notebook's discussion identifies the vertical impact (`θ = 0`) as the most demanding case. The plotted sweep supports the model's observed symmetric trend, with lower required stiffness for larger absolute impact angles.

> **Optimization caveat:** Bisection assumes the residual is appropriately bracketed over the stiffness interval; this is not checked explicitly in the notebook. The reported worst-case design is based on a finite angle grid rather than a proof over every continuous angle.

## Repository contents

| File | Description |
| --- | --- |
| [`Assignments_1.ipynb`](Assignments_1.ipynb) | Full simulation, POD, ROM validation, visualizations, and stiffness optimization |
| [`Assignments 1 - Problems description.pdf`](Assignments%201%20-%20Problems%20description.pdf) | Original assignment/problem description |

The notebook includes plots of the finite element displacement field, deformed domain, singular-value decay, POD modes, FOM–ROM comparison, and stiffness requirements against impact angle.

## Requirements and execution

The notebook is written in **Python** and uses:

- **FEniCS (legacy)** for finite element assembly [3]
- **`dlroms`** for finite element helper utilities (as imported in the notebook) [6]
- **NumPy**, **SciPy**, and **Matplotlib**
- **`ufl_legacy`** and a compatible FEniCS/UFL environment
- **Jupyter Notebook** or JupyterLab

The first cell attempts to import `dlroms` and, if unavailable, install it directly from its [GitHub repository](https://github.com/NicolaRFranco/dlroms). **This does not install FEniCS and all of its compiled dependencies.** Set up a working *legacy FEniCS* environment first; package compatibility depends on the chosen Python/FEniCS distribution. No pinned environment or dependency lockfile is provided in this repository.

Once the environment is ready, clone the repository, open [`Assignments_1.ipynb`](Assignments_1.ipynb), and run its cells in order.

```bash
jupyter notebook Assignments_1.ipynb
```

The notebook generates the snapshots and constructs the reduced model when executed; it does not load a pre-trained ROM from disk.

## Notes and limitations

- This is a **static linear-elastic approximation** of a hailstone load, not a time-dependent collision simulation.
- The ROM was evaluated against **50 held-out randomly sampled parameter combinations within the stated parameter ranges**; its accuracy outside this range was not established.
- The **~167×** timing comparison excludes the offline construction cost and may vary with hardware, software environment, and repeated measurements.
- The stiffness-selection sweep depends on the displacement threshold and load model used in the assignment; it should not be interpreted as a physically validated shield design.

## Techniques demonstrated

`Finite Element Method (FEM)` · `Linear Elasticity` · `Proper Orthogonal Decomposition (POD)` · `Singular Value Decomposition (SVD)` · `Reduced Basis Methods` · `Affine Decomposition` · `Offline–Online Model Reduction` · `Numerical Optimization` · `Python` · `FEniCS`


## References and citations

The **project-specific setup, code, numerical results, and optimization conclusions** come from [1] and [2]. References [3]–[6] credit the underlying scientific methods and software; they are **not** sources for the numerical results reported above.

1. **[1] Project implementation and results.** *Assignments_1.ipynb*, this repository. [View notebook](Assignments_1.ipynb).
2. **[2] Original problem statement.** *Assignments 1 - Problems description.pdf*, this repository. [View assignment](Assignments%201%20-%20Problems%20description.pdf).
3. **[3] FEniCS.** Alnæs, M. S., Blechta, J., Hake, J., et al. (2015). "The FEniCS Project Version 1.5." *Archive of Numerical Software*, **3**(100), 9–23. [https://doi.org/10.11588/ans.2015.100.20553](https://doi.org/10.11588/ans.2015.100.20553). See also the [official legacy FEniCS citation guide](https://fenicsproject.org/citing/legacy/).
4. **[4] Proper Orthogonal Decomposition.** Sirovich, L. (1987). "Turbulence and the Dynamics of Coherent Structures. Part I: Coherent Structures." *Quarterly of Applied Mathematics*, **45**(3), 561–571. [https://doi.org/10.1090/S0033-569X-1987-0910462-6](https://doi.org/10.1090/S0033-569X-1987-0910462-6).
5. **[5] Reduced-basis and affine decomposition.** Rozza, G., Huynh, D. B. P., & Patera, A. T. (2008). "Reduced Basis Approximation and a Posteriori Error Estimation for Affinely Parametrized Elliptic Coercive Partial Differential Equations: Application to Transport and Continuum Mechanics." *Archives of Computational Methods in Engineering*, **15**, 229–275. [https://doi.org/10.1007/s11831-008-9019-9](https://doi.org/10.1007/s11831-008-9019-9).
6. **[6] DLROMs software.** Franco, N. R. (2024). *NicolaRFranco/dlroms: First release (v1.0.0)*. Zenodo. [https://doi.org/10.5281/zenodo.13254758](https://doi.org/10.5281/zenodo.13254758). [Software repository](https://github.com/NicolaRFranco/dlroms).
