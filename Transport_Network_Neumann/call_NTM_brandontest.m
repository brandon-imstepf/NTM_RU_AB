% call_NTM_brandontest  Debug harness for the full RU-NTM + Amyloid-Beta
% coupling pipeline.
%
% Unlike call_NTM.m / call_NTM2.m, this script exposes every Aβ coupling
% parameter (r, Dif, kappa_mu_r, kappa_mu_u, kappa_delta, kappa_epsilon,
% coupling_term) in addition to the original kinetic/boundary-flux
% parameters, so the whole pipeline can be exercised end-to-end while
% debugging without hand-building a long name-value argument list each time.
%
% This is a SCRIPT, not a function: run it directly (F5 / "Run") and inspect
% N, M, R_Edge_Mass in the workspace afterward, or set breakpoints inside
% network_transport_model_Neu_func / NetworkFluxCalculator_Neu and step
% through with this as the entry point.
%
% All defaults below are chosen for FAST debug iteration (short time range,
% coarse mesh, small connectome subset) — NOT for a production/research run.
%
% NOTE on time range: network_transport_model_Neu_func's default 'trange'
% is NOT empty, so passing only 'dt'/'T' would silently be ignored (the
% code only falls back to dt/T when trange is empty). This script passes an
% explicit short trange directly to guarantee a fast run.

clear; clc;

% =========================================================================
% Data directory
% =========================================================================
matdir = 'C:\Users\USER\Documents\Code\NTM_ru_ab\Transport_Network_Neumann\MatFiles';

% =========================================================================
% Baseline kinetic / boundary-flux parameters (same meaning as call_NTM.m)
% =========================================================================
alpha_in    = 0;         % linear tau growth / recruitment rate
F_edge_0_in = 2;         % baseline distributed axonal source
beta_in     = 1e-06;     % aggregation saturation constant
gamma_in    = 0.001;     % aggregation rate gamma1
delta_in    = 50;        % motor-crowding factor (baseline, pre-AB-coupling)
epsilon_in  = 25;        % aggregation-velocity coupling (baseline, pre-AB-coupling)
lambda_in   = 1e-02;     % AIS diffusivity modifier

up_in_0  = 5;            % mu_u_0: uptake rate at soma end (baseline)
up_in_L  = 5;            % mu_u_L: uptake rate at synapse end (baseline)
rel_in_0 = 5;            % mu_r_0: release rate at soma end (baseline)
rel_in_L = 5;            % mu_r_L: release rate at synapse end (baseline)

init_rescale = 2e-2;     % target total tau (N+M) at t=0 in seeded regions
frac         = 0.92;     % free-tau fraction
axon_div     = 'r1';     % axon mass attribution

% =========================================================================
% Amyloid-Beta parameters
% =========================================================================
r_in   = 1;               % logistic growth rate for Abeta
Dif_in = 1;                % Abeta diffusion coefficient over the connectome

% --- Coupling selector: exactly ONE of these per run ---
%   'mu_u' | 'mu_r' | 'delta' | 'epsilon' | 'none'
coupling_term_in = 'mu_u';

kappa_mu_r_in    = 1;
kappa_mu_u_in    = 1;
kappa_delta_in   = 1;
kappa_epsilon_in = 1;

% =========================================================================
% Simulation setup — kept SMALL/FAST for debugging
% =========================================================================
connectome_subset_in = 'Hippocampus+PC+RSP';  % small region subset (~30 regions)
resmesh_in           = 'coarse';              % 250-pt axon mesh (vs 1000 for 'fine')
study_in             = 'Hurtado';

% Explicit short time range: 6 points (nt=6) => 5 time steps, all h<=50,
% so both the midpoint-predictor call and the per-step call get exercised.
trange_in = 0:0.01:0.05;

% =========================================================================
% Run the model
% =========================================================================
[N, M, R_Edge_Mass] = network_transport_model_Neu_func(matdir, ...
                                'alpha',             alpha_in,             ...
                                'F_edge_0',          F_edge_0_in,          ...
                                'beta',              beta_in,              ...
                                'gamma1',            gamma_in,             ...
                                'delta',             delta_in,             ...
                                'epsilon',           epsilon_in,           ...
                                'lambda1',           lambda_in,            ...
                                'mu_r_0',            rel_in_0,             ...
                                'mu_r_L',            rel_in_L,             ...
                                'mu_u_0',            up_in_0,              ...
                                'mu_u_L',            up_in_L,              ...
                                'init_rescale',      init_rescale,         ...
                                'frac',              frac,                 ...
                                'axon_div',          axon_div,             ...
                                'r',                 r_in,                 ...
                                'Dif',               Dif_in,               ...
                                'coupling_term',     coupling_term_in,     ...
                                'kappa_mu_r',        kappa_mu_r_in,        ...
                                'kappa_mu_u',        kappa_mu_u_in,        ...
                                'kappa_delta',       kappa_delta_in,       ...
                                'kappa_epsilon',     kappa_epsilon_in,     ...
                                'connectome_subset', connectome_subset_in, ...
                                'resmesh',           resmesh_in,           ...
                                'trange',            trange_in,            ...
                                'study',             study_in);

% =========================================================================
% Persist results, plus the coupling settings used, so each output file is
% self-describing (mirrors call_NTM2.m's convention of saving inputs too).
% =========================================================================
outdir       = fileparts(mfilename('fullpath'));
filename_out = fullfile(outdir, 'call_NTM_brandontest_output.mat');
save(filename_out, "N", "M", "R_Edge_Mass", "coupling_term_in", ...
     "kappa_mu_r_in", "kappa_mu_u_in", "kappa_delta_in", "kappa_epsilon_in");

fprintf('Saved debug run to %s\n', filename_out);
