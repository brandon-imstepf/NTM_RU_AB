function [network_flux_0, network_flux_L, mass_edge, F_source_edge, n_ss0, res_max, region_edge_mass] = ...
         NetworkFluxCalculator_Neu(tau_x0, tau_xL, f_ss0, Adj_input, varargin)
% NetworkFluxCalculator_Neu  Compute steady-state axonal fluxes across ALL network edges.
%
% =========================================================================
% OVERVIEW
% =========================================================================
% At each time step of the NTM, we assume that the intracellular tau
% distribution WITHIN every axon is at quasi-static steady state (the
% axonal PDE equilibrates much faster than the somatic ODE).  Given the
% current somatic tau concentrations (tau_x0, tau_xL) at both ends of
% each axon, this function:
%
%   1. For every connected pair (j→i) in the adjacency matrix, solves
%      a 1-D boundary-value problem for n(x) = intracellular tau density
%      along the axon.
%
%   2. Evaluates the resulting tau fluxes at x=0 (soma) and x=L (synapse).
%
%   3. Returns fluxes in matrices indexed as (source_region, target_region).
%
% The edge BVP is solved via a SHOOTING METHOD:
%   - The unknown B = n(0) (intracellular tau at the somatic end) is found
%     by fsolve so that the Neumann mass-balance condition at x=L is met.
%   - The per-edge solves are independent → parallelized with parfor.
%
% =========================================================================
% INPUTS
% =========================================================================
%   tau_x0     – soluble tau at the SOURCE (soma) end of each edge.
%                For the initial-time call: [1x3] array with representative
%                values for three edge-type classes.
%                For subsequent time steps: [nroi x nroi] matrix where
%                tau_x0(j,i) is tau at the soma of region j on the (j→i) axon.
%
%   tau_xL     – soluble tau at the TARGET (synapse) end of each edge.
%                For the initial-time call: [nroi x 1] vector (one per region).
%                For subsequent calls: [nroi x 1] vector (N at each region).
%
%   f_ss0      – initial guess for the shooting parameter n(0) for each edge.
%                Same shape as tau_x0 (used to warm-start fsolve for speed).
%
%   Adj_input  – full adjacency matrix from the CSV file [nroi_full x nroi_full].
%                Will be sub-set based on connectome_subset parameter.
%
%   varargin   – name-value pairs (see inputParser below)
%
% =========================================================================
% OUTPUTS
% =========================================================================
%   network_flux_0   – net tau flux at x=0 for each edge [nroi x nroi].
%                      network_flux_0(j,i) = flux on axon j→i at the soma of j.
%                      J_0 = -mu_r_0*n(0) + mu_u_0*tau_x0
%                      Positive = net tau departure from region j toward axon j→i.
%
%   network_flux_L   – net tau flux at x=L for each edge [nroi x nroi].
%                      network_flux_L(j,i) = flux on axon j→i at the synapse onto i.
%                      J_L = mu_r_L*n(L) - mu_u_L*tau_xL
%                      Positive = net tau arrival at region i from axon j→i.
%
%   mass_edge        – total intracellular tau mass on each axon (currently 0;
%                      computation is commented out for performance).
%
%   F_source_edge    – cumulative distributed source F_edge_0*(x3-x0) on each axon [nroi x nroi].
%                      Represents tau produced along the axon itself.
%
%   n_ss0            – solved shooting parameter n(0) for each edge [nroi x nroi].
%                      Passed back as the warm-start guess for the next time step.
%
%   res_max          – maximum absolute residual of the shooting equation across
%                      all edges (quality-of-solve indicator; should be ~0).
%
%   region_edge_mass – breakdown of edge mass attributed to source vs. target region
%                      (currently 0; placeholder for regional mass analysis).

% =========================================================================
% SECTION 1: Default parameter values
% (same biophysics as the standalone solvers; see those for full derivations)
% =========================================================================
beta_         = 1e-06;
gamma1_x0_    = 2e-05;    % aggregation rate at source end
gamma1_xL_    = 2e-05;    % aggregation rate at target end
gamma2_       = 0;
delta_        = 50;
delta_x0_     = 0.01;
delta_xL_     = 0.01;
epsilon_      = 25;
epsilon_x0_   = 0.025;
epsilon_xL_   = 0.025;
lambda1_x0_   = 0.0025;   % AIS permeability at source end
lambda1_xL_   = 0.0025;   % AIS permeability at target end
frac_         = 0.92;     % diffusing fraction of tau (Konsack 2007)
L_int_        = 1000;     % axon length in um
L1_           = 200;      % somatodendritic length in um
L_ais_        = 40;       % AIS length in um
L_syn_        = 40;       % synaptic terminal length in um
resmesh_      = 'coarse'; % 'fine' or 'coarse' spatial mesh
reltol_       = 1e-6;
abstol_       = 1e-6;
fsolvetol_    = 1e-6;
connectome_subset_ = 'Hippocampus+PC+RSP';
len_scale_    = 1e-3;     % um → mm
time_scale_   = 1;        % set to 1 here; caller may pass actual time scale

% AB-Specific Parameters (Brandon, 7/31/2026)
% To do


% Boundary flux constants
F_edge_0_      = 1e-08 * time_scale_;
mu_r_0_        = 1e-05 * time_scale_;   % release rate at soma end
mu_r_L_        = mu_r_0_;               % release rate at synapse end
mu_u_0_        = 1e-04 * time_scale_;   % uptake rate at soma end
mu_u_L_        = mu_u_0_;               % uptake rate at synapse end
F_edge_dt      = 0e-05;
F_edge_dx      = 0e-05;
idx_netw_x0_   = 1;   % default: edge is active at x=0
idx_netw_xL_   = 1;   % default: edge is active at x=L
axon_div_      = 'r1';

% =========================================================================
% SECTION 2: inputParser
% =========================================================================
ip = inputParser;
validScalar  =      @(x) isnumeric(x) && isscalar(x) && (x >= 0);
validArray   =      @(x) isnumeric(x);
validLogical =      @(x) islogical(x);
validAxonDiv =      @(x) strcmp(x,'r1') | strcmp(x,'half') | strcmp(x,'r2');


addParameter(ip, 'beta',              beta_,          validScalar);
addParameter(ip, 'F_edge_0',          F_edge_0_,      validScalar);
addParameter(ip, 'gamma2',            gamma2_,        validScalar);
addParameter(ip, 'delta',             delta_,         validScalar);
addParameter(ip, 'delta_x0',          delta_x0_,      validArray);
addParameter(ip, 'delta_xL',          delta_xL_,      validArray); % Added Arrays for Delta, Epsilon (Brandon, 8/10/2026)
addParameter(ip, 'epsilon',           epsilon_,       validScalar);
addParameter(ip, 'epsilon_x0',        epsilon_x0_,    validArray);
addParameter(ip, 'epsilon_xL',        epsilon_xL_,    validArray);
addParameter(ip, 'frac',              frac_,          validScalar);
addParameter(ip, 'gamma1_x0',         gamma1_x0_,     validArray);
addParameter(ip, 'gamma1_xL',         gamma1_xL_,     validArray);
addParameter(ip, 'lambda1_x0',        lambda1_x0_,    validArray);
addParameter(ip, 'lambda1_xL',        lambda1_xL_,    validArray);
addParameter(ip, 'idx_netw_x0',       idx_netw_x0_,   validArray);
addParameter(ip, 'idx_netw_xL',       idx_netw_xL_,   validLogical);
addParameter(ip, 'mu_r_0',            mu_r_0_,        validArray); % Changed mu's to Array for diffuse support (Brandon, 7/31/2026)
addParameter(ip, 'mu_r_L',            mu_r_L_,        validArray); 
addParameter(ip, 'mu_u_0',            mu_u_0_,        validArray);
addParameter(ip, 'mu_u_L',            mu_u_L_,        validArray);
addParameter(ip, 'L_int',             L_int_,         validScalar);
addParameter(ip, 'L1',                L1_,            validScalar);
addParameter(ip, 'resmesh',           resmesh_);
addParameter(ip, 'L_ais',             L_ais_);
addParameter(ip, 'L_syn',             L_syn_);
addParameter(ip, 'reltol',            reltol_,        validScalar);
addParameter(ip, 'abstol',            abstol_,        validScalar);
addParameter(ip, 'fsolvetol',         fsolvetol_,     validScalar);
addParameter(ip, 'connectome_subset', connectome_subset_);
addParameter(ip, 'len_scale',         len_scale_,     validScalar);
addParameter(ip, 'time_scale',        time_scale_,    validScalar);
addParameter(ip, 'axon_div',          axon_div_,      validAxonDiv);
parse(ip, varargin{:});

% =========================================================================
% SECTION 3: Derived physical constants (scaled to model units)
% =========================================================================
beta_new   = ip.Results.beta * ip.Results.time_scale;

% Convert lengths from um to model units
L1_new    = ip.Results.L1    * ip.Results.len_scale;
L_int_new = ip.Results.L_int * ip.Results.len_scale;
L_ais_new = ip.Results.L_ais * ip.Results.len_scale;
L_syn_new = ip.Results.L_syn * ip.Results.len_scale;

% Biophysical velocities and diffusivity (Konsack 2007, scaled)
v_a    = 0.7 * ip.Results.len_scale * ip.Results.time_scale;   % anterograde velocity
v_r    = 0.7 * ip.Results.len_scale * ip.Results.time_scale;   % retrograde velocity
diff_n = 12  * ip.Results.len_scale^2 * ip.Results.time_scale; % free tau diffusivity

% Boundary flux rates as Arrays
mu_r_0 = ip.Results.mu_r_0;
mu_r_L = ip.Results.mu_r_L;
mu_u_0 = ip.Results.mu_u_0;
mu_u_L = ip.Results.mu_u_L;

% Extract gamma1 and lambda1 arrays (can be per-region scalars or matrices)
gamma1_x0  = ip.Results.gamma1_x0;
gamma1_xL  = ip.Results.gamma1_xL;
lambda1_x0 = ip.Results.lambda1_x0;
lambda1_xL = ip.Results.lambda1_xL;
idx_netw_x0 = ip.Results.idx_netw_x0;
idx_netw_xL = ip.Results.idx_netw_xL;

% Extract delta, epsilon arrays
delta_x0 = ip.Results.delta_x0;
delta_xL = ip.Results.delta_xL;
epsilon_x0 = ip.Results.epsilon_x0;
epsilon_xL = ip.Results.epsilon_xL;

% =========================================================================
% SECTION 4: Build the spatial mesh
%
% The axon domain is [0, L1+L_int-L_syn] (presyn + AIS + axon).
% Finer spacing is used near the presyn/AIS boundary (x = L1) and near
% the AIS/axon boundary (x = L1+L_ais) to capture rapid changes in n(x).
%
% Two resolution modes:
%   'fine'   – 1000 mesh points (used for final results)
%   'coarse' – 250 mesh points  (used for faster parameter sweeps)
% =========================================================================
if strcmp(ip.Results.resmesh, 'fine')
    num_comp = 1000;
    num_ext  = 100;    % points in somatodendritic region
    num_int  = num_comp - 2 * num_ext;
    xmesh1 = [linspace(0, L1_new - 10*ip.Results.len_scale, num_ext-40), ...
               (L1_new - 9.75*ip.Results.len_scale) : 0.25*ip.Results.len_scale : L1_new];
    xmesh_int = [(L1_new + 0.25*ip.Results.len_scale) : 0.25*ip.Results.len_scale : (L1_new + L_ais_new), ...
                  linspace(L1_new + L_ais_new + 0.25*ip.Results.len_scale, L1_new + L_int_new - L_syn_new, ...
                           num_int - ((L_ais_new + L_syn_new) / (0.25*ip.Results.len_scale))), ...
                 (L1_new + L_int_new - (L_syn_new - 0.25*ip.Results.len_scale)) : 0.25*ip.Results.len_scale : (L1_new + L_int_new - L_syn_new)];
elseif strcmp(ip.Results.resmesh, 'coarse')
    num_comp = 250;
    num_ext  = 25;
    num_int  = num_comp - 2 * num_ext;
    xmesh1 = [linspace(0, L1_new - 10*ip.Results.len_scale, num_ext-5), ...
               (L1_new - 8*ip.Results.len_scale) : 2*ip.Results.len_scale : L1_new];
    xmesh_int = [(L1_new + 2*ip.Results.len_scale) : 2*ip.Results.len_scale : (L1_new + L_ais_new), ...
                  linspace(L1_new + L_ais_new + 2*ip.Results.len_scale, L1_new + L_int_new - L_syn_new, ...
                           num_int - ((L_ais_new + L_syn_new) / (2*ip.Results.len_scale))), ...
                 (L1_new + L_int_new - (L_syn_new - 2*ip.Results.len_scale)) : 2*ip.Results.len_scale : (L1_new + L_int_new - L_syn_new)];
end
xmesh = [xmesh1, xmesh_int];

% =========================================================================
% SECTION 5: Set up the compartment-specific steady-state ODEs
%
% All three ODEs are constructed as lambda functions (closures) that capture
% the current parameters.  They will be evaluated inside the parfor loop.
%
% Note: these are defined at the outer scope so that the parfor body can
% call them.  MATLAB closures are copy-captured in parfor, which is fine.
% =========================================================================

% Averaging functions for gamma1, lambda1, delta & epsilon along an edge (i→j):
%   Use the arithmetic mean of the source-end and target-end values.
gamma1_fun   = @(gamma1_i, gamma1_j, x) (gamma1_i + gamma1_j) ./ 2;
lambda_fun   = @(lambda_i, lambda_j, x) (lambda_i + lambda_j) ./ 2;
delta_fun    = @(delta_i, delta_j, x) (delta_i + delta_j) ./ 2;
epsilon_fun  = @(epsilon_i, epsilon_j, x) (epsilon_i + epsilon_j) ./ 2;

% Distributed source: int_F_edge(x) = F_edge_0 * (x - x0)
x0 = xmesh(1);
int_F_edge = @(x, idx) (ip.Results.F_edge_0 .* (x - x0)) .* idx;

% --- Spatial masks: logical index vectors into xmesh ---
presyn_mask = (xmesh <= L1_new);
xmesh_presyn = xmesh(presyn_mask);

ais_mask = logical(-1 + (xmesh > L1_new) + (xmesh < (L1_new + L_ais_new)));
xmesh_ais = xmesh(ais_mask);

axon_mask = logical(-1 + (xmesh >= L1_new + L_ais_new) + ...
                    (xmesh <= (L1_new + L_int_new - L_syn_new)));
xmesh_axon = xmesh(axon_mask);

% Compartment boundary points (used as ICs for chained ODEs)
x1 = xmesh_presyn(end);  % presyn / AIS boundary
x2 = xmesh_ais(end);     % AIS / axon boundary
x3 = xmesh_axon(end);    % synaptic end of axon

% ODE solver options
options = odeset('RelTol', ip.Results.reltol, 'AbsTol', ip.Results.abstol);

% --- Compartment 1: Somatodendritic region [0, L1] ---
% dn/dx = (1/D)*(mu_r_0*B - mu_u_0*N_i - int_F_edge(x))
% B = n(0) = shooting parameter; N_i = extracellular tau at source node
n0 = @(B) B;   % IC: n(x=0) = B
n_ss_presyn = @(B, N_i, idx,mu_r_0_i, mu_u_0_i) ode15s(@(x, n) ode_ss_n(x, B, n, diff_n, N_i, mu_r_0_i, mu_u_0_i, idx), ...
                                      [0, L1_new], n0(B), options);
n_ss_presyn = @(B, N_i, idx,mu_r_0_i, mu_u_0_i, x) deval(n_ss_presyn(B, N_i, idx,mu_r_0_i, mu_u_0_i), x);

% --- Compartment 2: AIS [L1, L1+L_ais] ---
% Same ODE as presyn but diffusivity scaled by lambda (AIS permeability factor):
%   D_ais = diff_n * lambda_fun(lambda_i, lambda_j)
% IC: continuity at x1 (n must be equal across presyn/AIS boundary)
n_ss_ais = @(B, N_i, idx, lambda1_i, lambda1_j, mu_r_0_i, mu_u_0_i) ...
    ode15s(@(x, n) ode_ss_n(x, B, n, diff_n * lambda_fun(lambda1_i, lambda1_j, x), N_i, mu_r_0_i, mu_u_0_i, idx), ...
           [L1_new, L1_new + L_ais_new], n_ss_presyn(B, N_i, idx, mu_r_0_i, mu_u_0_i, x1), options);
n_ss_ais = @(B, N_i, idx, lambda1_i, lambda1_j, mu_r_0_i, mu_u_0_i, x) deval(n_ss_ais(B, N_i, idx, lambda1_i, lambda1_j,mu_r_0_i, mu_u_0_i), x);

% --- Compartment 3: Axon proper [L1+L_ais, L1+L_int-L_syn] ---
% Nonlinear ODE including active transport:
%   dn/dx = (1/(frac*D)) * [(1-frac)*v_eff(n)*n + mu_r_0*B - mu_u_0*N_i - F_edge(x)]
% where v_eff captures motor crowding (delta) and aggregation coupling (epsilon).
% IC: continuity at x2 (end of AIS)
n_ss_axon = @(B, N_i, idx, lambda1_i, lambda1_j, gamma1_i, gamma1_j, delta_i,delta_j,epsilon_i,epsilon_j,mu_r_0_i, mu_u_0_i) ...
    ode15s(@(x, n) ode_ss_axon(x, B, gamma1_i, gamma1_j, delta_i, delta_j, epsilon_i, epsilon_j, n, N_i, mu_r_0_i, mu_u_0_i, idx), ...
           [L1_new + L_ais_new, L1_new + L_int_new - L_syn_new], ...
           n_ss_ais(B, N_i, idx, lambda1_i, lambda1_j, mu_r_0_i, mu_u_0_i, x2), options);
n_ss_axon = @(B, N_i, idx, lambda1_i, lambda1_j, gamma1_i, gamma1_j, delta_i, delta_j, epsilon_i, epsilon_j, mu_r_0_i, mu_u_0_i, x) ...
    deval(n_ss_axon(B, N_i, idx, lambda1_i, lambda1_j, gamma1_i, gamma1_j, delta_i, delta_j, epsilon_i, epsilon_j,mu_r_0_i, mu_u_0_i), x);

% --- Shooting residual f_init(B) ---
% Mass-balance condition (Neumann BC at x=L):
%   -mu_r_0*B + mu_u_0*N_i + int_F_edge(x3) - mu_r_L*n(x3) + mu_u_L*N_j = 0
% Meaning: flux in at x=0 + distributed production = flux out at x=L
f_init = @(B, N_i, N_j, idx, lambda1_i, lambda1_j, gamma1_i, gamma1_j, delta_i, delta_j, epsilon_i, epsilon_j,mu_r_0_i, mu_u_0_i,mu_r_L_i, mu_u_L_i) ...
    -mu_r_0_i .* B + mu_u_0_i .* N_i + int_F_edge(x3, idx) ...
    - mu_r_L_i .* n_ss_axon(B, N_i, idx, lambda1_i, lambda1_j, gamma1_i, gamma1_j, delta_i, delta_j, epsilon_i, epsilon_j, mu_r_0_i, mu_u_0_i, x3) ...
    + mu_u_L_i .* N_j;

% =========================================================================
% SECTION 6: Sub-set the adjacency matrix to the selected brain region group
%
% The full mouse atlas has 426 regions, but we may only model a subset
% for speed (e.g. hippocampus + piriform + retrosplenial cortex = ~30 regions).
% =========================================================================
Adj = Adj_input;

switch ip.Results.connectome_subset
    case 'Hippocampus'
        Adj = Adj([27:37 (27+213):(37+213)], [27:37 (27+213):(37+213)]);
    case 'Hippocampus+PC+RSP'
        % Hippocampus (rows 27-37), Piriform (78-80), RSP (147), and LH mirrors
        adjinds = [27:37, 78:80, 147];
        adjinds = [adjinds, adjinds + 213];
        Adj = Adj(adjinds, adjinds);
    case 'RH'
        Adj = Adj(1:213, 1:213);    % right hemisphere only
    case 'LH'
        Adj = Adj(214:end, 214:end); % left hemisphere only
    case 'Single'
        Adj = 1;                     % single-edge test
    case 'Single_bis'
        Adj = [1 1 1];              % three-edge test (seed↔seed, seed↔non-seed)
end
nroi = size(Adj, 2);  % number of regions in the chosen subset

% Normalize mu_r/mu_u to full Adj-shaped arrays if the caller passed scalars
% (driver still passes flat scalars today (8/27/26); will pass Mu_r_hat/Mu_u_hat arrays later)
if isscalar(mu_r_0), mu_r_0 = mu_r_0 * ones(size(Adj)); end
if isscalar(mu_u_0), mu_u_0 = mu_u_0 * ones(size(Adj)); end
if isscalar(mu_r_L), mu_r_L = mu_r_L * ones(nroi, 1);   end
if isscalar(mu_u_L), mu_u_L = mu_u_L * ones(nroi, 1);   end
% debugging 8/27/2026
assert(isequal(size(mu_r_0), size(Adj)), 'mu_r_0 shape does not match Adj'); % rwk_debug
assert(isequal(size(mu_u_0), size(Adj)), 'mu_u_0 shape does not match Adj'); % rwk_debug
assert(numel(mu_r_L) == nroi, 'mu_r_L length does not match nroi'); % rwk_debug
assert(numel(mu_u_L) == nroi, 'mu_u_L length does not match nroi'); % rwk_debug

% =========================================================================
% SECTION 7: Allocate output arrays
% All matrices are [nroi x nroi]: entry (j,i) corresponds to the edge j→i.
% Only entries where Adj(j,i) = 1 will be filled; the rest remain 0.
% =========================================================================
network_flux_0 = zeros(size(Adj));
network_flux_L = zeros(size(Adj));
n_ss0          = zeros(size(Adj));
res            = zeros(size(Adj));
F_source_edge  = zeros(size(Adj));

% Temporary "dict" arrays used inside parfor (parfor cannot write to sliced
% shared arrays directly when the slice indices depend on the loop variable;
% the dict approach avoids the issue by collecting all outputs first, then
% unpacking them in a serial loop afterward).
n_ss0_dict          = zeros(size(Adj));
res_dict            = zeros(size(Adj));
F_source_edge_dict  = zeros(size(Adj));
network_flux_0_dict = zeros(size(Adj));
network_flux_L_dict = zeros(size(Adj));

if isscalar(mu_r_0), mu_r_0 = mu_r_0 * ones(size(Adj)); end
if isscalar(mu_u_0), mu_u_0 = mu_u_0 * ones(size(Adj)); end
if isscalar(mu_r_L), mu_r_L = mu_r_L * ones(nroi, 1);   end
if isscalar(mu_u_L), mu_u_L = mu_u_L * ones(nroi, 1);   end
if isscalar(delta_x0),   delta_x0   = delta_x0   * ones(size(Adj)); end
if isscalar(epsilon_x0), epsilon_x0 = epsilon_x0 * ones(size(Adj)); end
if isscalar(delta_xL),   delta_xL   = delta_xL   * ones(nroi, 1);   end
if isscalar(epsilon_xL), epsilon_xL = epsilon_xL * ones(nroi, 1);   end
% debugging 8/27/2026
assert(isequal(size(mu_r_0), size(Adj)), 'mu_r_0 shape does not match Adj'); % rwk_debug
assert(isequal(size(mu_u_0), size(Adj)), 'mu_u_0 shape does not match Adj'); % rwk_debug
assert(numel(mu_r_L) == nroi, 'mu_r_L length does not match nroi'); % rwk_debug
assert(numel(mu_u_L) == nroi, 'mu_u_L length does not match nroi'); % rwk_debug
assert(isequal(size(delta_x0), size(Adj)), 'delta_x0 shape does not match Adj'); % rwk_debug
assert(isequal(size(epsilon_x0), size(Adj)), 'epsilon_x0 shape does not match Adj'); % rwk_debug
assert(numel(delta_xL) == nroi, 'delta_xL length does not match nroi'); % rwk_debug
assert(numel(epsilon_xL) == nroi, 'epsilon_xL length does not match nroi'); % rwk_debug

% =========================================================================
% SECTION 8: Parallel loop over target regions i
%
% For each target region i, we solve the steady-state axon BVP for EVERY
% source region j that has a connection to i (i.e., Adj(j,i) = 1).
%
% This loop is parallelized via parfor (requires Parallel Computing Toolbox).
% Each iteration is independent: the steady-state n(x) on axon j→i depends
% only on tau_x0(j,i) and tau_xL(i), not on any other edge.
%
% The shooting solve is warm-started from f_ss0(j,i) (the n(0) found at
% the previous time step) so fsolve converges quickly.
%
% parfor restriction: Adj, n_ss0_dict, etc. must be indexed by i as full
% rows/columns (sliced variables), not by arbitrary index sets, so we
% accumulate into temp arrays and then unpack below.
% =========================================================================
parfor i = 1:nroi

    % ------------------------------------------------------------------
    % Identify which source regions j are connected to target i.
    % Adj_in is a logical column vector: Adj_in(j) = true if j→i exists.
    % ------------------------------------------------------------------
    Adj_in = logical(Adj(:, i));  % [nroi x 1] logical

    % Replicate the scalar target-end quantities to match the number of sources
    tau_xL_i       = repmat(tau_xL(i),       length(Adj_in), 1);
    idx_netw_xL_i  = repmat(idx_netw_xL(i),  length(Adj_in), 1);
    gamma1_xL_i    = repmat(gamma1_xL(i),    length(Adj_in), 1);
    lambda1_xL_i   = repmat(lambda1_xL(i),   length(Adj_in), 1);
    delta_xL_i     = repmat(delta_xL(i),     length(Adj_in), 1);
    epsilon_xL_i   = repmat(epsilon_xL(i),   length(Adj_in), 1);
    mu_r_L_i       = repmat(mu_r_L(i),       length(Adj_in), 1);
    mu_u_L_i       = repmat(mu_u_L(i),       length(Adj_in), 1);

    % Apply the adjacency mask: keep only entries where j is connected to i
    i_app = Adj_in;  % indices of active source regions for this target i

    tau_x0_i      = tau_x0(i_app, i);      % source-end tau for active edges
    tau_xL_i      = tau_xL_i(i_app);        % target-end tau (scalar broadcast)
    idx_netw_x0_i = idx_netw_x0(i_app);     % is source seeded?
    idx_netw_xL_i = idx_netw_xL_i(i_app);   % is target seeded?
    gamma1_x0_i   = gamma1_x0(i_app, i);   % aggregation rate at source
    gamma1_xL_i   = gamma1_xL_i(i_app);     % aggregation rate at target
    lambda1_x0_i  = lambda1_x0(i_app, i);  % AIS permeability at source
    lambda1_xL_i  = lambda1_xL_i(i_app);    % AIS permeability at target
    delta_x0_i    = delta_x0(i_app, i);
    delta_xL_i    = delta_xL_i(i_app);
    epsilon_x0_i  = epsilon_x0(i_app,i);
    epsilon_xL_i  = epsilon_xL_i(i_app);
    mu_r_L_i      = mu_r_L_i(i_app);
    mu_u_L_i      = mu_u_L_i(i_app);
    mu_r_0_i      = mu_r_0(i_app, i);
    mu_u_0_i      = mu_u_0(i_app, i);

    assert(numel(mu_r_0_i) == sum(i_app), 'mu_r_0_i length mismatch vs i_app'); % rwk_debug
    assert(numel(mu_r_L_i) == sum(i_app), 'mu_r_L_i length mismatch vs i_app'); % rwk_debug

    if ~isempty(tau_x0_i)  % skip if region i has no incoming connections
        options_fsolve = optimset('TolFun', ip.Results.fsolvetol, 'Display', 'off');

        % idx_netw: edge is "active" if EITHER endpoint has seeded tau.
        % This controls whether F_edge (distributed source) is applied.
        idx_netw = logical(idx_netw_x0_i + idx_netw_xL_i);

        % Build the shooting residual for all source regions simultaneously
        fun_ss = @(B) f_init(B, tau_x0_i, tau_xL_i, idx_netw, ...
                              lambda1_x0_i, lambda1_xL_i, gamma1_x0_i, gamma1_xL_i, ...
                              delta_x0_i, delta_xL_i, epsilon_x0_i, epsilon_xL_i,mu_r_0_i,mu_u_0_i,mu_r_L_i,mu_u_L_i);

        % Warm-start fsolve from the previous step's n(0) values (+ tiny perturbation)
        f0 = f_ss0(i_app, i) + 1e-7;

        % Solve for the shooting parameter n(0) on all active (j→i) axons
        n_ss0_temp = fsolve(fun_ss, f0, options_fsolve);

        % ------------------------------------------------------------------
        % Pack results into full-nroi arrays (zero for inactive edges).
        % We write to temp arrays here because parfor cannot directly update
        % slices of n_ss0_dict(i,:) where i_app selects non-contiguous rows.
        % ------------------------------------------------------------------
        temp_n_ss0_arr = zeros(nroi, 1);
        temp_n_ss0_arr(i_app) = n_ss0_temp;
        n_ss0_dict(i, :) = temp_n_ss0_arr;

        % Shooting residual (should be ~machine epsilon if converged)
        res_temp = fun_ss(n_ss0_temp);
        temp_res_arr = zeros(nroi, 1);
        temp_res_arr(i_app) = res_temp;
        res_dict(i, :) = temp_res_arr;

        % Distributed source on each active edge:
        %   F_source_edge = int_F_edge(x3) = F_edge_0 * (x3 - x0)
        F_source_edge_temp = int_F_edge(x3, idx_netw);
        temp_F_source = zeros(nroi, 1);
        temp_F_source(i_app) = F_source_edge_temp;
        F_source_edge_dict(i, :) = temp_F_source;

        % ------------------------------------------------------------------
        % Flux at x=0 (soma end of source region j):
        %   J_0 = -mu_r_0 * n(0) + mu_u_0 * N_source
        % Positive J_0 → net tau flow from source node INTO the axon.
        % ------------------------------------------------------------------
        netflux_0_temp = -mu_r_0_i .* n_ss0_temp + mu_u_0_i .* tau_x0_i;
        temp_netflux_0 = zeros(nroi, 1);
        temp_netflux_0(i_app) = netflux_0_temp;
        network_flux_0_dict(i, :) = temp_netflux_0;

        % ------------------------------------------------------------------
        % Flux at x=L (synapse end → arriving at target region i):
        %   J_L = mu_r_L * n(x3) - mu_u_L * N_target
        % Positive J_L → net tau delivery from axon INTO target node.
        % n(x3) is evaluated by re-integrating the full axon ODE with
        % the just-solved n(0) = n_ss0_temp.
        % ------------------------------------------------------------------
        netflux_L_temp = mu_r_L_i .* n_ss_axon(n_ss0_temp, tau_x0_i, idx_netw, ...
                                               lambda1_x0_i, lambda1_xL_i, ...
                                               gamma1_x0_i, gamma1_xL_i, ...
                                               delta_x0_i, delta_xL_i, ...
                                               epsilon_x0_i, epsilon_xL_i, ...
                                               mu_r_0_i, mu_u_0_i, x3) ...
                         - mu_u_L_i .* tau_xL_i;
        temp_netflux_L = zeros(nroi, 1);
        temp_netflux_L(i_app) = netflux_L_temp;
        network_flux_L_dict(i, :) = temp_netflux_L;
    end
end  % parfor over target regions i

% =========================================================================
% SECTION 9: Unpack the parfor results into the output matrices
%
% This serial loop copies from the row-indexed "dict" arrays into the
% column-indexed output matrices.  The asymmetry arises because parfor
% iterated over i (target), but the output matrices are indexed as (j,i)
% (source j, target i) and we wrote n_ss0_dict(i,j) = value for edge j→i.
% =========================================================================
for i = 1:nroi
    Adj_in = logical(Adj(:, i));
    i_app  = Adj_in;

    n_ss0(i_app, i)         = n_ss0_dict(i, i_app);
    res(i_app, i)           = res_dict(i, i_app);
    F_source_edge(i_app, i) = F_source_edge_dict(i, i_app);
    network_flux_0(i_app, i) = network_flux_0_dict(i, i_app);
    network_flux_L(i_app, i) = network_flux_L_dict(i, i_app);
end

% =========================================================================
% SECTION 10: Edge mass (currently disabled; placeholder)
%
% The intracellular tau mass on each axon (∫n dx + ∫m dx) would be
% computed here and attributed to source and/or target region based on
% axon_div.  Commented out for performance; R_Edge_Mass is returned as 0
% from the top-level model.
% =========================================================================
mass_edge        = 0;
region_edge_mass = 0;

% Maximum residual across all edges: a quality-of-convergence indicator.
% Should be near machine epsilon if fsolve converged everywhere.
res_max = max(abs(res), [], 'all');


% =========================================================================
% NESTED FUNCTIONS
% These are defined at the outer function scope and captured by the closures
% n_ss_presyn, n_ss_ais, n_ss_axon defined above.
% =========================================================================

    % ODE for somatodendritic and AIS compartments (pure linear diffusion)
    % State: n = intracellular mobile tau density
    % RHS:  dn/dx = (1/D) * (release - uptake - distributed_production)
    %   mu_r_0*B  : release at soma (proportional to n(0)=B; constant per solve)
    %   mu_u_0*N_i: uptake of extracellular N_i into the axon
    %   int_F_edge: cumulative distributed source term
    function nprime = ode_ss_n(x, B, n, D, N_i, mu_r_0_i, mu_u_0_i, idx)
        nprime = (1./D) .* (mu_r_0_i .* B - mu_u_0_i .* N_i - int_F_edge(x, idx));
    end

    % ODE for the axon proper (nonlinear: includes active transport)
    % State: n = intracellular mobile tau density
    % RHS includes diffusive term + advective (active-transport) contribution:
    %   dn/dx = (1/(frac*D)) * [(1-frac)*v_eff(n)*n + source terms]
    % where v_eff = v_a*(1+delta*n)*(1 - gamma1*epsilon*n^2/(beta-gamma2*n)) - v_r
    %   v_a*(1+delta*n)   : anterograde velocity enhanced by motor crowding
    %   *(1-...)           : speed reduction when aggregated tau (∝ n^2) clogs motors
    %   -v_r              : net retrograde component subtracted
    function nprime = ode_ss_axon(x, B, gamma1_i, gamma1_j, delta_i, delta_j, epsilon_i, epsilon_j, n, N_i, mu_r_0_i, mu_u_0_i, idx)
        nprime = 1 ./ (ip.Results.frac * diff_n) .* ...
                 ((1 - ip.Results.frac) .* ...
                  ((v_a .* (1 + delta_fun(delta_i,delta_j,x) .* n) .* ...
                    (1 - (gamma1_fun(gamma1_i, gamma1_j,x) .* epsilon_fun(epsilon_i,epsilon_j,x) .* n.^2) ./ ...
                         (beta_new - ip.Results.gamma2 .* n)) - v_r)) .* n + ...
                  mu_r_0_i .* B - mu_u_0_i .* N_i - int_F_edge(x, idx));
        % debugging 8/27/2026
        d_avg = delta_fun(delta_i, delta_j,x);
        e_avg = epsilon_fun(epsilon_i, epsilon_j,x);
        assert(isequal(size(d_avg),size(n)) || isscalar(d_avg), 'ode_ss_axon: delta avg shape mismatch'); % rwk_debug
        assert(isequal(size(e_avg),size(n)) || isscalar(e_avg), 'ode_ss_axon: epsilon avg shape mismatch'); % rwk_debug
    end

end  % NetworkFluxCalculator_Neu