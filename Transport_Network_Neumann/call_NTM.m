function [] = call_NTM(matdir, filename_in, alpha_in, F_edge_0_in, beta_in, ...
                        gamma_in, delta_in, epsilon_in, lambda_in, up_in_0, up_in_L, ...
                        rel_in_0, rel_in_L, init_rescale, frac, axon_div)
% call_NTM  Entry-point wrapper: run the Network Transport Model and save results.
%
% This function is a thin shell around network_transport_model_Neu_func.
% Its only job is to accept named parameter values, forward them to the
% model, and write the output to a .mat file.
%
% The underlying model simulates how misfolded tau protein spreads through
% the mouse brain connectome over time.  Tau is represented by two
% populations at each brain-region node:
%   N  – soluble (extracellular/interstitial) tau  [nroi x nt]
%   M  – bound/aggregated tau, in quasi-static equilibrium with N  [nroi x nt]
%
% INPUTS
%   matdir        – path to the folder containing the required .mat data
%                   files (connectome, atlas, mouse tauopathy data, etc.)
%   filename_in   – full path of the output .mat file to create/overwrite
%
%   --- Biological / kinetic parameters ---
%   alpha_in      – linear growth rate for tau at each node (recruitment).
%                   Units: 1/time_scale.  Set to 0 to disable.
%   F_edge_0_in   – baseline axonal production rate (distributed source of
%                   tau along every active axon edge). Units: concentration/length/time.
%   beta_in       – saturation constant in the tau-aggregation equilibrium
%                   M = gamma1*N^2 / (beta - gamma2*N).  Prevents M from
%                   diverging; represents the finite capacity of the
%                   aggregation pathway.
%   gamma_in      – aggregation rate gamma1: controls how fast soluble tau
%                   converts into the bound pool M.
%   delta_in      – motor-crowding enhancement factor.  In the axon ODE,
%                   the anterograde velocity scales as v_a*(1 + delta*n),
%                   meaning higher intracellular tau density speeds up
%                   anterograde transport (saturating motors are recruited).
%   epsilon_in    – aggregation-velocity coupling.  High epsilon means that
%                   aggregated tau strongly slows anterograde transport
%                   (the term -(gamma1*epsilon*n^2)/(beta-gamma2*n) in v_eff).
%   lambda_in     – AIS (axon initial segment) diffusivity modifier.
%                   Diffusivity in the AIS = frac*D*lambda, so lambda < 1
%                   creates a diffusion barrier at the AIS.
%
%   --- Boundary flux parameters (at each axon endpoint) ---
%   up_in_0       – mu_u_0: uptake rate of soluble tau N_i into the axon
%                   at the somatic end (x = 0).
%   up_in_L       – mu_u_L: uptake rate of soluble tau N_j into the axon
%                   at the synaptic end (x = L).
%   rel_in_0      – mu_r_0: release rate of intracellular axonal tau n(0)
%                   back into the extracellular space at the somatic end.
%   rel_in_L      – mu_r_L: release rate of intracellular axonal tau n(L)
%                   into the extracellular space at the synaptic end.
%
%   --- Simulation setup ---
%   init_rescale  – target value of N + M at t=0 in seeded regions.
%                   Used to find the initial N via fsolve (since M depends
%                   non-linearly on N through the aggregation formula).
%   frac          – fraction of intracellular tau that is in the *diffusing*
%                   (free) state vs. actively transported.  From Konsack 2007
%                   the default is 0.92.
%   axon_div      – string specifying how axon mass is attributed to source
%                   vs. target region: 'r1' (all to source), 'r2' (all to
%                   target), or 'half' (split at midpoint).
%
% OUTPUTS  (saved to filename_in)
%   N             – soluble tau concentration, size [nroi x nt]
%   M             – bound tau concentration, size [nroi x nt]
%   Reg_edge_mass – regional edge mass (currently zeroed out; placeholder)

% -------------------------------------------------------------------------
% Call the main model function.
% All parameters are forwarded as name-value pairs so network_transport_model_Neu_func
% can validate them via inputParser and apply defaults for any omitted ones.
% -------------------------------------------------------------------------
[N, M, Reg_edge_mass] = network_transport_model_Neu_func(matdir, ...
                                          'alpha',        alpha_in,    ...
                                          'F_edge_0',     F_edge_0_in, ...
                                          'beta',         beta_in,     ...
                                          'gamma1',       gamma_in,    ...
                                          'delta',        delta_in,    ...
                                          'epsilon',      epsilon_in,  ...
                                          'lambda1',      lambda_in,   ...
                                          'mu_r_0',       rel_in_0,    ...
                                          'mu_r_L',       rel_in_L,    ...
                                          'mu_u_0',       up_in_0,     ...
                                          'mu_u_L',       up_in_L,     ...
                                          'init_rescale', init_rescale, ...
                                          'frac',         frac,        ...
                                          'axon_div',     axon_div);

% -------------------------------------------------------------------------
% Persist results.  N and M are the primary outputs used for downstream
% analysis (e.g., comparison with mouse tauopathy data).
% -------------------------------------------------------------------------
save(filename_in, "N", "M", "Reg_edge_mass");
