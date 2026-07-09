% run_script.m  Minimal development/testing script.
%
% Calls network_transport_model_Neu_func with all default parameter values
% (defined inside that function) and saves the output to a hard-coded path.
%
% Use this script to do a quick end-to-end test of the model without
% having to specify every parameter.  For production runs or parameter
% sweeps, use call_NTM or NTM_sim instead.

% -------------------------------------------------------------------------
% Run the model with default parameters.
% No matdir is passed, so network_transport_model_Neu_func will look for
% the MatFiles folder inside the current working directory.
% -------------------------------------------------------------------------
[N, M] = network_transport_model_Neu_func();

%% Save variables
% N = soluble tau [nroi x nt]
% M = bound/aggregated tau [nroi x nt]
save("C:/Users/USER/Documents/Code/NTM_ru/output/savefile.mat", "N", "M");
