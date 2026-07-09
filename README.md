# Network Transport Model (NTM) — Tau Propagation in the Mouse Brain

## Overview

This codebase implements the **Network Transport Model (NTM)**, a mathematical model of how misfolded tau protein spreads through the mouse brain connectome. It is relevant to Alzheimer's disease and other tauopathies, where tau aggregates progressively invade connected brain regions via axonal transport.

The model operates at two coupled spatial scales:

1. **Brain-network scale** — each brain region is a node; tau concentration evolves in time via a system of ODEs.
2. **Single-axon scale** — each axonal connection between two regions is resolved as a 1-D PDE, solved at steady state to compute how much tau flows along that axon.

---

## Directory Structure

```
NTM_ru_ab/
├── Transport_Network_Neumann/     # All MATLAB source files
│   ├── MatFiles/                  # Required data files (not in repo)
│   ├── call_NTM.m
│   ├── call_NTM2.m
│   ├── NTM_sim.m
│   ├── run_script.m
│   ├── network_transport_model_Neu_func.m
│   ├── NetworkFluxCalculator_Neu.m
│   ├── DataToCCF.m
│   ├── Flux_Calculator_RU.m
│   └── Flux_Calculator_Neu_edge.m
├── output/                        # Simulation output .mat files
└── README.md
```

### Required Data Files (`MatFiles/`)

| File | Contents |
|---|---|
| `Mouse_Tauopathy_Data_HigherQ.mat` | Measured tau data from mouse experiments (e.g. Hurtado study) |
| `DefaultAtlas.mat` | Brain region volumes from the Allen CCF atlas |
| `CCF_labels.mat` | Region names, parent regions, hemisphere labels (426 regions) |
| `Connectomes.mat` | Axonal connectivity matrices |
| `mouse_adj_matrix_19_01.csv` | Binary adjacency matrix for the full connectome |

---

## Call Chain

```
NTM_sim.m  ──────────────────────────────────────┐
run_script.m  ────────────────────────────────────┤
call_NTM.m  ──────────────────────────────────────┤──► network_transport_model_Neu_func.m
call_NTM2.m  ─────────────────────────────────────┘         │
                                                             ├──► NetworkFluxCalculator_Neu.m
                                                             └──► DataToCCF.m

Standalone diagnostics (not called during a real run):
    Flux_Calculator_RU.m
    Flux_Calculator_Neu_edge.m
```

---

## File Descriptions

### Entry Points

#### `call_NTM.m`
The primary entry point for production runs. Accepts all model parameters as function arguments, calls `network_transport_model_Neu_func`, and saves `N`, `M`, and `Reg_edge_mass` to a `.mat` file.

#### `call_NTM2.m`
A simplified entry point used by `NTM_sim` for parameter sweeps. Hard-codes several parameters (beta, lambda, frac, etc.) and exposes only the parameters being varied. Also saves the input parameter values alongside the outputs for traceability.

#### `NTM_sim.m`
A Monte Carlo / parameter sweep harness. Given an integer random seed, it draws one set of model parameters from uniform distributions over biologically plausible ranges, then calls `call_NTM2`. Designed to be called in parallel (e.g. via a SLURM array job) with different seeds.

| Parameter | Range |
|---|---|
| `alpha` (growth rate) | [0, 1] |
| `gamma1` (aggregation rate) | [0.001, 0.008] |
| `delta` (motor crowding factor) | [10, 100] |
| `epsilon` (aggregation-velocity coupling) | [10, 100] |
| `uprel` (combined uptake/release rate) | [50, 778] |

#### `run_script.m`
A minimal two-line script that runs the model with all default parameters. Useful for a quick end-to-end test.

---

### Core Model

#### `network_transport_model_Neu_func.m`
The main solver. Takes all model parameters as name-value pairs, runs the full simulation, and returns `N` and `M`.

**What it does, step by step:**

1. **Loads data** — connectome, atlas volumes, CCF labels, mouse tauopathy data.
2. **Preprocesses the connectome** — zeros the diagonal, thresholds weak connections at 80% of the mean non-zero weight.
3. **Sets up the initial condition** — maps experimental seed data onto CCF regions via `DataToCCF`, then computes the initial soluble tau `N(0)` by solving `N + M(N) = init_rescale` with `fsolve`.
4. **Subsets the connectome** — restricts to a chosen brain-region group (e.g. Hippocampus + Piriform + RSP, or the full graph).
5. **Computes initial fluxes** — calls `NetworkFluxCalculator_Neu` with three representative edge types to get `J_0` and `J_L` at `t=0`, then broadcasts these onto the full flux matrices.
6. **Time-steps** — loops over `t(1)…t(nt)`, updating `N` at each step.
7. **Computes `M`** — derives bound tau from the quasi-static equilibrium formula.

**Time integration scheme:**

| Steps | Method | Reason |
|---|---|---|
| h ≤ 50 | Midpoint / RK2 | Higher accuracy during the rapidly-changing early phase |
| h > 50 | Forward Euler | Cheaper; acceptable once N is varying smoothly |

---

### Axon Flux Engine

#### `NetworkFluxCalculator_Neu.m`
Called at every time step to compute the steady-state axonal flux on every active edge in the network. This is the most computationally expensive part of the model.

**What it does:**

For each target region `i`, it loops (in parallel via `parfor`) over all source regions `j` connected to `i`, and for each active edge `j→i`:

1. Solves the 1-D axon BVP (three compartments) using a shooting method.
2. Returns the flux at `x=0` (soma) and `x=L` (synapse).

**Outputs** (all `[nroi × nroi]` matrices, indexed as `(source, target)`):
- `network_flux_0` — flux at soma end of each axon
- `network_flux_L` — flux at synapse end of each axon
- `n_ss0` — shooting parameter `n(0)` for each edge (used as warm-start at next step)
- `F_source_edge` — cumulative distributed tau source on each axon
- `res_max` — maximum shooting residual (convergence quality indicator)

---

### Data Utility

#### `DataToCCF.m`
Maps experimental tau data from study-specific ROIs (typically ~50–100 regions) onto the 426-region Allen CCF atlas. A single study ROI can correspond to several CCF sub-regions (many-to-one mapping). CCF regions not covered by the study are set to `NaN`.

---

### Standalone Diagnostics

These two files are **not called during a real simulation run**. They exist for development, parameter exploration, and generating plots of the axon PDE solution.

#### `Flux_Calculator_RU.m`
An earlier version of the single-axon edge solver. Uses a full flux-at-L Neumann boundary condition formula. Takes a quasi-static time `time_qs` as its first argument.

#### `Flux_Calculator_Neu_edge.m`
A newer single-axon diagnostic that uses the same simplified "global mass balance" Neumann BC as the production solver. Generates a plot of `n(x)` along the axon.

---

## Mathematical Background

### Two Tau Populations

At each brain region node `i`:

- `N_i(t)` — soluble (extracellular/interstitial) tau. This is the tau that can be taken up by axons and transported.
- `M_i(t)` — aggregated (bound) tau. Assumed to be in quasi-static equilibrium with `N`:

```
M_i = γ₁ · N_i² / (β - γ₂ · N_i)
```

### Network ODE

The full time derivative is:

```
d/dt [N_i + M_i] = (1/Vol_i) · [Σ_j Conn_ji · J_L(j→i)  −  Σ_j Conn_ij · J_0(i→j)]  +  α · N_i
```

Using the chain rule `dM/dt = (dM/dN) · dN/dt`:

```
(1 + m_t,i) · dN_i/dt = (1/Vol_i) · [flux_in_i − flux_out_i]  +  α · N_i
```

where `m_t,i = dM/dN = γ₁ · N_i · (2β − γ₂·N_i) / (β − γ₂·N_i)²`

### Single-Axon BVP (Spatial Compartments)

Each axon `j→i` has three spatial compartments along `x ∈ [0, L]`:

```
[0,       L₁]             Somatodendritic (SD) region  — pure diffusion
[L₁,      L₁+L_ais]       Axon initial segment (AIS)   — diffusion, reduced by λ₁
[L₁+L_ais, L₁+L_int-L_syn] Axon proper                 — diffusion + active transport
```

The steady-state ODE in each compartment:

**SD and AIS** (linear diffusion):
```
dn/dx = (1/D) · (μ_r0 · n(0) − μ_u0 · N_source − F_edge(x))
```

**Axon proper** (nonlinear — includes active transport):
```
dn/dx = (1/(f·D)) · [(1−f) · v_eff(n) · n  +  μ_r0 · n(0) − μ_u0 · N_source − F_edge(x)]
```

where the effective transport velocity is:
```
v_eff(n) = v_a · (1 + δ·n) · (1 − γ₁·ε·n² / (β − γ₂·n))  −  v_r
```
- `v_a·(1 + δ·n)` — anterograde velocity enhanced by motor crowding (`δ`)
- The factor `(1 − ...)` — speed reduction when aggregated tau clogs motor proteins (`ε`)
- `−v_r` — net retrograde component subtracted

### Boundary Conditions and Fluxes

The unknown `n(0) = B` (intracellular tau at the soma end) is found by `fsolve` to satisfy the **global mass-balance condition**:

```
−μ_r0 · B  +  μ_u0 · N_source  +  F_edge_total  −  μ_rL · n(L)  +  μ_uL · N_target  =  0
```

The resulting boundary fluxes are:

```
J_0 = −μ_r0 · n(0) + μ_u0 · N_source      (positive → tau leaves source node into axon)
J_L =  μ_rL · n(L) − μ_uL · N_target      (positive → tau arrives at target node from axon)
```

---

## Key Parameters

| Parameter | Symbol | Typical Value | Meaning |
|---|---|---|---|
| `alpha` | α | 0 | Linear tau growth / recruitment rate |
| `beta` | β | 1×10⁻⁶ | Aggregation saturation constant |
| `gamma1` | γ₁ | 0.001 | Aggregation rate (soluble → bound) |
| `delta` | δ | 50 | Motor crowding enhancement factor |
| `epsilon` | ε | 25 | Aggregation-velocity coupling |
| `lambda1` | λ₁ | 0.01 | AIS diffusivity modifier (barrier strength) |
| `frac` | f | 0.92 | Fraction of tau in the diffusing (free) pool |
| `F_edge_0` | — | 2 | Baseline distributed tau source along axons |
| `mu_r_0` | μ_r0 | 5 | Release rate of intracellular tau at soma |
| `mu_u_0` | μ_u0 | 5 | Uptake rate of extracellular tau at soma |
| `init_rescale` | — | 0.02 | Target total tau (N+M) in seeded regions at t=0 |
| `L_int` | — | 1000 μm | Axon length |
| `L1` | — | 200 μm | Somatodendritic region length |
| `L_ais` | — | 40 μm | Axon initial segment length |

---

## Outputs

Both `N` and `M` are `[nroi × nt]` matrices:

- **Rows** — brain regions (in the order of `CCF_labels` after subsetting)
- **Columns** — time points (from the `trange` vector, in model time units of months)

`N(:, end)` gives the final soluble tau distribution; compare against mouse tauopathy data to evaluate model fit.
