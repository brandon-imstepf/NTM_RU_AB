function [] = NTM_sim(matdir, seed_in, study_in)
% NTM_sim  Monte-Carlo / parameter-sweep harness for the Network Transport Model.
%
% Draws a set of model parameters from uniform distributions (one draw per
% call, seeded by seed_in for reproducibility) and then runs call_NTM2 to
% simulate tau propagation through the mouse connectome.
%
% Intended use: call this function many times with different seed_in values
% (e.g. from a SLURM array job) to build a distribution of model outputs
% across the biologically plausible parameter space.
%
% INPUTS
%   matdir   – path to folder containing the .mat data files required by
%              the model (connectome, atlas, tauopathy data).
%   seed_in  – integer random seed.  Controls the parameter draw so each
%              call is deterministic given the same seed.
%   study_in – string selecting which mouse tauopathy study seeds tau at t=0
%              (e.g. 'Hurtado').  Passed through to call_NTM2.
%
% OUTPUT
%   Writes a .mat file named  simulation_<seed_in>.mat  to the hard-coded
%   output path below.  Each file contains N, M, and the parameter values
%   used for that run (see call_NTM2 for the full list).

% -------------------------------------------------------------------------
% NOTE: matdir is overridden here with a hard-coded path.
% If you move the data files, update this line.
% -------------------------------------------------------------------------
matdir = 'C:\Users\USER\Documents\Code\NTM_ru\Transport_Network_Neumann\MatFiles';

% -------------------------------------------------------------------------
% Seed MATLAB's random number generator so this draw is reproducible.
% Using the same seed_in will always produce the same parameter set.
% -------------------------------------------------------------------------
rng(seed_in);

% -------------------------------------------------------------------------
% Draw one uniform random number per free parameter.
% Each r_* is in [0, 1] and is linearly mapped to a biologically plausible range.
% -------------------------------------------------------------------------
r_alpha   = rand();   % maps to alpha   range [0, 1]
r_gamma   = rand();   % maps to gamma1  range [0.001, 0.008]
r_delta   = rand();   % maps to delta   range [10, 100]
r_epsilon = rand();   % maps to epsilon range [10, 100]
r_uprel   = rand();   % maps to mu_r/u  range [50, 778]

% -------------------------------------------------------------------------
% Scale the draws into their biologically motivated parameter ranges.
% Formula: value = (max - min)*r + min
% -------------------------------------------------------------------------
alpha_in   = r_alpha;                            % growth rate [0, 1]
gamma_in   = ((0.008 - 0.001) * r_gamma) + 0.001; % aggregation rate gamma1 [0.001, 0.008]
delta_in   = ((100   - 10)    * r_delta) + 10;    % motor crowding factor [10, 100]
epsilon_in = ((100   - 10)    * r_epsilon) + 10;  % aggregation-velocity coupling [10, 100]
uprel_in   = ((778   - 50)    * r_uprel) + 50;    % combined uptake/release rate [50, 778]

% -------------------------------------------------------------------------
% Build the output file path.  Each simulation gets a unique file named by
% its seed so runs can be submitted in parallel without file conflicts.
% -------------------------------------------------------------------------
save_file = "/Users/nbarron/Desktop/simulation_" + num2str(seed_in) + ".mat";

% -------------------------------------------------------------------------
% Run the model with the sampled parameters.
% Results (N, M, plus parameter values) are written to save_file.
% -------------------------------------------------------------------------
call_NTM2(matdir, save_file, study_in, alpha_in, gamma_in, delta_in, epsilon_in, uprel_in);
