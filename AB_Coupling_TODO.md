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

**Debugging convention:** code new sections with debugging in mind —
liberal `fprintf` statements to trace intermediate values, and MATLAB
`assert(cond, 'message')` calls to check invariants (e.g. array shapes
line up, values stay in expected ranges) as we go, not just after something
breaks. Brandon appends his own trailing comment tag to each debug
line/block so they're greppable and safe to strip once a section is
confirmed working (e.g. searching the file for that tag before a "clean"
commit). Claude should proactively suggest specific `fprintf`/`assert`
additions alongside any new code guidance from here on, sized to what's
actually risky in that step (shape mismatches, unexpected NaNs/Infs,
sign errors) rather than blanket-instrumenting everything.

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
- [x] **`Mu_u_hat`** built: `ip.Results.mu_u_0 + ip.Results.kappa_mu_u * Abeta` (line 434)
- [x] **`Mu_r_hat`** built: `ip.Results.mu_r_0 + ip.Results.kappa_mu_r * Abeta` (line 437)
- [x] **`Delta_hat`** built: `ip.Results.delta + ip.Results.kappa_delta * Abeta` (line 440)
- [x] **`Epsilon_hat`** built: `ip.Results.epsilon + ip.Results.kappa_epsilon * Abeta` (line 443)
- [x] **`gamma1`/aggregation: resolved — leaving `Gamma1` as-is, no Aβ coupling
      for now.** Replacing it with a `Gamma1_hat` built from `Abeta` would have
      required reworking `Gamma_app` (the midpoint approximation used in the
      h≤50 RK2 branch), since `Gamma_app` only works today because the old
      ROI-ramp has a closed form evaluable at any time — `Abeta` is only
      tabulated at the `nt` grid points from its `ode45` solve, with no
      closed-form midpoint value. Existing `gamma_time_dip_roi` ramp mechanism
      stays untouched and Aβ-independent.
- [x] **Resolved — `delta`'s sign:** `delta` multiplies the anterograde (`v_a`)
      term, not `v_r`. Proposal Table 1 calls it "Retrograde enhancement
      factor" — inconsistent with the actual mechanism; user will raise with
      advisor separately. Coupling itself stays additive/positive, no sign flip.
- [x] **`NetworkFluxCalculator_Neu.m` array-plumbing bugs — all fixed** (were
      open at last check-in, confirmed fixed on this pass):
      - `repmat` syntax error at the `delta_xL_i`/`epsilon_xL_i` lines (~389-390) — fixed.
      - `epsilon_x0`/`epsilon_xL` validators changed `validScalar` → `validArray` — fixed.
      - Missing semicolon on `epsilon_fun` (~248) — fixed.
- [x] **Delta/epsilon threading — done.** Full chain verified end-to-end
      (`fun_ss` → `f_init` → `n_ss_axon` → `ode_ss_axon`, all argument counts
      and orders cross-checked at every hop): `delta_i/delta_j/epsilon_i/
      epsilon_j` flow from the `parfor` slicing through to `v_eff` via
      `delta_fun`/`epsilon_fun` in place of the old bare `ip.Results.delta`/
      `ip.Results.epsilon`.
- [x] **`mu_r_0`/`mu_u_0`/`mu_r_L`/`mu_u_L` now per-edge sliced — done.**
      Added a shape-normalization guard (~line 338-348, right after `nroi`/
      `Adj` are known) that broadcasts a scalar caller input to
      `size(Adj)`/`[nroi 1]` so the function stays backward-compatible with
      today's driver (which still passes flat scalars) while also supporting
      real per-region arrays later. `mu_r_0_i`/`mu_u_0_i` (via `(i_app,i)`)
      and `mu_r_L_i`/`mu_u_L_i` (via `repmat`+mask, mirroring `gamma1_xL_i`)
      threaded through `n_ss_presyn`→`n_ss_ais`→`n_ss_axon`→`ode_ss_n`/
      `ode_ss_axon`, `f_init`, `fun_ss`, and the direct `netflux_0_temp`/
      `netflux_L_temp` computations. Several argument-order bugs surfaced and
      were fixed during review (idx/mu cross-wiring in `ode_ss_n`'s two
      callers, `n`/`N_i` vs `mu_r_0_i`/`mu_u_0_i` swap in `ode_ss_axon`'s
      call, missing args in `n_ss_ais`'s call into `n_ss_presyn`, missing
      delta/epsilon args in the `netflux_L_temp` call into `n_ss_axon`) —
      all confirmed fixed on final pass. `rwk_debug`-tagged asserts added at
      the shape-normalization guard, the `parfor` slicing, and inside
      `ode_ss_axon`. **Not yet run** — next step is an actual test call to
      see whether the asserts pass and results look physically sane.
- [x] **`coupling_term` selector — done, correctly placed in the driver.**
      Initially added by mistake to `NetworkFluxCalculator_Neu.m` (which has
      no `Abeta`/`nt`/`kappa_*` in scope at all — it's a per-timestep, per-edge
      BVP solver, architecturally agnostic to Aβ); moved to
      `network_transport_model_Neu_func.m`'s `inputParser`
      (`'mu_u'|'mu_r'|'delta'|'epsilon'|'none'`, validated via
      `validCouplingTerm`; fixed a missing-comma bug in that validator —
      `strcmp(x {...})` parsed as cell-indexing, not two args). Section 10's
      `Mu_u_hat`/`Mu_r_hat`/`Delta_hat`/`Epsilon_hat` builds are now guarded
      by it, `else` branch fills a full `[nroi x nt]` constant so downstream
      shape-dependent indexing keeps working regardless of which term is hot.
- [x] **All three `NetworkFluxCalculator_Neu` call sites wired — done.**
      - Init call (t=0): intentionally left passing flat `ip.Results.*`
        scalars — consistent with its own pre-existing `gamma1_x0=gamma1_new`
        simplification (uses representative edge types, not real region
        identities, so there's no meaningful per-region hat value to plug in).
      - Per-step call (every iteration, after the h≤50/h>50 branch resolves
        `N(:,h+1)`): passes `Mu_r_hat`/`Mu_u_hat`/`Delta_hat`/`Epsilon_hat`
        sliced at `h+1` (`Mu_r_0_h1 = Mu_r_hat(:,h+1).*Adj` etc., mirroring
        `Gamma1_x0_h1`). Fixed a missing-paren syntax error found here
        (`Mu_u_hat(:,h+1,` — swallowed the rest of the call's arguments).
      - Midpoint-predictor call (h≤50 branch only): same pattern but sliced
        at `h` (not `h+1`, since `N(:,h+1)` doesn't exist yet at the
        predictor stage) — `Mu_r_0_app = Mu_r_hat(:,h).*Adj` etc.
      **As of now, Aβ coupling actually reaches the flux solver end-to-end,**
      gated correctly by `coupling_term`.
- [x] **Debug harness added: `call_NTM_brandontest.m`** (Transport_Network_Neumann/).
      A script (not a function, unlike `call_NTM`/`call_NTM2`) exposing every
      Aβ parameter (`r`, `Dif`, `kappa_mu_r`, `kappa_mu_u`, `kappa_delta`,
      `kappa_epsilon`, `coupling_term`) alongside the original kinetic params.
      Runs a short (5-step) fast debug pass — small connectome subset, coarse
      mesh, explicit short `trange` (passing only `dt`/`T` would silently be
      ignored, since the model's default `trange` is non-empty and takes
      priority) — hitting all three call sites and all `rwk_debug` asserts.
      Saves to `call_NTM_brandontest_output.mat`. **Not yet actually run** —
      this is the next step: confirm the whole pipeline executes cleanly and
      the asserts pass, then try each `coupling_term` value in turn to confirm
      the selector actually isolates one term at a time.

---

## Open Questions To Revisit
1. `gamma2` (secondary/seeded aggregation) is 0 throughout — confirmed intentional,
   keeps `M(N)` invertible for the ADNI `N_i(0)` recovery formula. No action needed,
   just don't reintroduce it without revisiting that formula.
2. `Abeta` is currently spatially uniform (varies only in time), because
   `B0 = 0.05*ones(nroi,1)` is uniform and the graph-Laplacian diffusion term
   is identically zero when all regions start (and stay) equal — flagged to
   user, no decision made yet on whether to seed `B0` non-uniformly instead.
