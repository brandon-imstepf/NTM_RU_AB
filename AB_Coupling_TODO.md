# RU-NTM + Amyloid-Beta Coupling — Project Tracker

## Overall Goal

Extend the Release-Uptake Network Transport Model (RU-NTM) in this repo with
**state-dependent Amyloid-Beta (Aβ) coupling**, per the PhD proposal
"Extending the Release–Uptake Network Transport Model with State-Dependent
Amyloid-β Coupling."

Aβ is modeled as a standalone field `B_i(t)` per brain region, evolving
independently of tau (one-way coupling: Aβ influences tau kinetics, tau does
not influence Aβ) via:

```
dB_i/dt = r*B_i*(1 - B_i) + Dif * sum_j[ C_ij*(B_j - B_i) ]
```

(logistic growth + Aβ diffusion over the connectome graph Laplacian). This is
pre-solved once via `ode45` before the tau time-stepping loop begins, and the
resulting `Abeta(i, h)` trajectory is then used to modulate tau kinetics
parameters at each node/time step. **Corrected convention (additive, not
multiplicative):** each "hat" quantity is the existing baseline parameter
plus a coupling term, built directly from the parameters already in the
model — not from the separately-declared `_hat` scalar parameters (those are
being left declared but unused):

- **Release**: `Mu_r_hat = mu_r_0 + kappa_mu_r*Abeta`
- **Uptake**: `Mu_u_hat = mu_u_0 + kappa_mu_u*Abeta`
- **Aggregation** (`gamma1`): **deferred — see Status below.** `Gamma1` is
  being left as-is (no Aβ coupling for now); avoids the `Gamma_app` midpoint
  problem (no closed form for `Abeta` between time grid points).
- **Retrograde/anterograde transport** (`delta`): `Delta_hat = delta + kappa_delta*Abeta`
  (base stays positive, "suppressing" route — see resolved note below)
- **Aggregation velocity** (`epsilon`): new addition, not yet scoped —
  `Epsilon_hat = epsilon + kappa_epsilon*Abeta` (parameter scaffolding for
  `kappa_epsilon` likely still needs to be added)

Longer-term (not started): Sobol sensitivity analysis over the new coupling
parameters, ADNI human data comparison (using the `N_i(0) = sqrt(beta*M_i(0)/(gamma1*B_i(0)))`
recovery formula), and identifiability analysis per proposal Section 3.2.

**Working rule for this project:** the user (Brandon) writes all code changes
himself. Claude's role is strictly review/debug/explain — never implement.

---

## Status

### Step 1 — Aβ standalone pre-solve — done
- [x] `ode45` pre-solve added in Section 8.5 of `network_transport_model_Neu_func.m`
- [x] Correctly placed after Section 7 (connectome subsetting) and Section 8 (time vector `t`)
- [x] `r` parameter added and threaded via `ip.Results.r`
- [x] **Fixed:** `Dif` now scales the whole Laplacian term
      (`ip.Results.Dif*(Conn*B - sum(Conn,2).*B)`), confirmed by user.
- [x] Coupling parameter pairs declared in `inputParser`. Note: the original
      `mu_r_hat`/`mu_u_hat`/`gamma_hat`/`delta_hat` scalar defaults (all = 1)
      are now **dead/unused** — superseded by the additive convention below,
      which builds `_hat` arrays straight from existing baseline params
      (`mu_r_0`, `mu_u_0`, `delta`, `gamma1_new`). User is leaving the dead
      scalars declared but untouched.

### Step 2 — Wire coupling terms into the actual flux/aggregation computations — in progress
- [x] **`Mu_u_hat`** built: `ip.Results.mu_u_0 + ip.Results.kappa_mu_u * Abeta` (line 432)
- [x] **`Mu_r_hat`** built: `ip.Results.mu_r_0 + ip.Results.kappa_mu_r * Abeta` (line 435,
      fixed — was accidentally multiplicative, now matches `Mu_u_hat`'s form)
- [ ] **`Mu_r_hat`/`Mu_u_hat` still need array support in `NetworkFluxCalculator_Neu.m`**
      before they can take effect — currently `mu_r_0/mu_r_L/mu_u_0/mu_u_L` are
      `validScalar`-only and captured as bare outer-scope scalars in nested
      closures (not threaded as explicit arguments, unlike `gamma1_x0/gamma1_xL`).
- [x] **`gamma1`/aggregation: resolved — leaving `Gamma1` as-is, no Aβ coupling
      for now.** Replacing it with a `Gamma1_hat` built from `Abeta` would have
      required reworking `Gamma_app` (the midpoint approximation used in the
      h≤50 RK2 branch), since `Gamma_app` only works today because the old
      ROI-ramp has a closed form evaluable at any time — `Abeta` is only
      tabulated at the `nt` grid points from its `ode45` solve, with no
      closed-form midpoint value. User decided this cost isn't worth it right
      now; the existing `gamma_time_dip_roi` ramp mechanism stays untouched.
      Revisit if a later step needs gamma1 coupled for consistency.
- [ ] **`delta`:** needs the most new plumbing — no array support exists at all
      today (`ip.Results.delta` used as a bare scalar in `v_eff`). Plan: add
      `delta_x0`/`delta_xL` array params (validated via `validArray`, same
      pattern as `gamma1_x0`/`gamma1_xL`), a `delta_fun` averaging helper
      mirroring `gamma1_fun = @(a,b,x) (a+b)./2`, and thread through `v_eff`/
      `ode_ss_axon` in place of the current bare scalar. `Delta_hat` array
      itself not yet built in the driver.
- [x] **Resolved — `delta`'s sign:** `delta` multiplies the anterograde (`v_a`)
      term, not `v_r` (retrograde has no tunable in this codebase at all).
      Proposal Table 1 calls `delta` "Retrograde enhancement factor," which is
      inconsistent with its actual mechanism — user plans to raise this
      modeling gap with their advisor separately. For the coupling itself:
      **decided to keep `Delta_hat` positive**, additive form
      `delta + kappa_delta*Abeta` — no sign flip, no new `v_r`-specific
      multiplier. Wiring plan unchanged: add `delta_x0`/`delta_xL` array
      params + `delta_fun` (mirroring `gamma1_fun`), thread through
      `v_eff`/`ode_ss_axon`.
- [ ] **`epsilon` (new scope addition):** user intends an `Epsilon_hat` as well
      (`epsilon + kappa_epsilon*Abeta`). No `kappa_epsilon` parameter exists
      yet — needs default + `addParameter` in Section 1/2. `epsilon` is
      currently `validScalar`-only in `NetworkFluxCalculator_Neu.m`, used raw
      in `v_eff`'s aggregation-quench term — same array-plumbing gap as `delta`.

---

## Open Questions To Revisit
1. `gamma2` (secondary/seeded aggregation) is 0 throughout — confirmed intentional,
   keeps `M(N)` invertible for the ADNI `N_i(0)` recovery formula. No action needed,
   just don't reintroduce it without revisiting that formula.
2. Once `Mu_r_hat`/`Mu_u_hat`/`Delta_hat`/`Epsilon_hat` are all built in the driver,
   `NetworkFluxCalculator_Neu.m` needs matching array-support additions for
   `mu_r_0/mu_r_L/mu_u_0/mu_u_L`, `delta`, and `epsilon` (currently all
   `validScalar`-only/bare-scalar) before any of these couplings actually affect
   the flux computation.
