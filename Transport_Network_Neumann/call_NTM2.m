function [] = call_NTM2(matdir, filename_in, study_in, alpha_in, gamma_in, ...
                         delta_in, epsilon_in, uprel_in)
% call_NTM2  Simplified entry-point wrapper used by NTM_sim for parameter sweeps.
%
% Like call_NTM, this is a thin shell around network_transport_model_Neu_func.
% Compared to call_NTM it:
%   - Exposes fewer free parameters (mu_r and mu_u are tied together as uprel_in)
%   - Hard-codes several less-varied parameters (beta, lambda1, frac, axon_div, etc.)
%   - Accepts a 'study_in' argument to select which mouse tauopathy experiment
%     provides the seed distribution and atlas data
%   - Saves additional metadata (the input parameter values) alongside N and M
%
% This function is typically called by NTM_sim, which samples the free
% parameters randomly before calling here.
%
% INPUTS
%   matdir      – path to folder with .mat data files (connectome, atlas, etc.)
%   filename_in – full path of the .mat file to write results into
%   study_in    – string key selecting a study from mousedata_struct
%                 (e.g. 'Hurtado', 'IbaP301S').  Determines which brain
%                 regions are seeded with tau at t=0.
%   alpha_in    – linear tau growth / recruitment rate (0 = no growth)
%   gamma_in    – aggregation rate gamma1 (soluble -> bound conversion)
%   delta_in    – anterograde motor crowding factor
%   epsilon_in  – aggregation-velocity coupling factor
%   uprel_in    – combined uptake/release rate at both axon endpoints.
%                 Applied as: mu_r_0 = mu_r_L = mu_u_0 = mu_u_L = uprel_in
%                 (i.e. symmetric release and uptake at soma and synapse).
%
% OUTPUTS  (saved to filename_in)
%   N           – soluble tau [nroi x nt]
%   M           – bound/aggregated tau [nroi x nt]
%   alpha_in, gamma_in, delta_in, epsilon_in, uprel_in, study_in
%               – the parameter values used, saved for bookkeeping

% -------------------------------------------------------------------------
% Hard-coded defaults for parameters that are not varied in this sweep.
% -------------------------------------------------------------------------
F_edge_0_in = 0;        % No distributed axonal tau production source
beta_in     = 0.000001; % Aggregation saturation constant
lambda_in   = 0.025;    % AIS diffusivity modifier (25% of axon diffusivity)
init_rescale = 0.02;    % Target total tau (N+M) in seeded regions at t=0
frac        = 0.7;      % Fraction of tau in diffusing state (vs. actively transported)
axon_div    = 'half';   % Split axon mass equally between source and target region
conn_subset = 'all';    % Use the full connectome (all brain regions)
% conn_subset = 'Hippocampus+PC+RSP';  % Alternative: hippocampus + piriform + RSP subset

% -------------------------------------------------------------------------
% Run the model.  uprel_in is used as ALL four boundary flux parameters so
% that the release and uptake rates are identical at both axon endpoints.
% -------------------------------------------------------------------------
[N, M, ~] = network_transport_model_Neu_func(matdir, ...
                                          'alpha',             alpha_in,    ...
                                          'F_edge_0',          F_edge_0_in, ...
                                          'beta',              beta_in,     ...
                                          'gamma1',            gamma_in,    ...
                                          'delta',             delta_in,    ...
                                          'epsilon',           epsilon_in,  ...
                                          'lambda1',           lambda_in,   ...
                                          'mu_r_0',            uprel_in,    ...
                                          'mu_r_L',            uprel_in,    ...
                                          'mu_u_0',            uprel_in,    ...
                                          'mu_u_L',            uprel_in,    ...
                                          'init_rescale',      init_rescale, ...
                                          'frac',              frac,        ...
                                          'axon_div',          axon_div,    ...
                                          'connectome_subset', conn_subset, ...
                                          'study',             study_in,    ...
                                          'resmesh',           'coarse');   % coarse mesh for speed

% -------------------------------------------------------------------------
% Save results together with the parameter values so each output file is
% self-describing and can be traced back to its inputs.
% -------------------------------------------------------------------------
save(filename_in, "N", "M", "alpha_in", "gamma_in", "delta_in", ...
                  "epsilon_in", "uprel_in", "study_in");
