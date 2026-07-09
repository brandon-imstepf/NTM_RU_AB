function [J_0, J_L] = Flux_Calculator_Neu_edge(varargin)
% Flux_Calculator_Neu_edge  Standalone diagnostic for a single axon edge.
%
% PURPOSE
%   A self-contained solver for the steady-state intracellular tau profile
%   along a single axon (edge).  Used for development, parameter exploration,
%   and generating Figure-quality plots of n(x).
%
% RELATIONSHIP TO OTHER FILES
%   This is a NEWER version of Flux_Calculator_RU that uses the "Neumann"
%   (release-uptake) boundary-condition formulation adopted throughout
%   the full network model (NetworkFluxCalculator_Neu).
%   The key simplification vs. Flux_Calculator_RU: the mass-balance BC at
%   x=L is written as:
%       -mu_r_0*B + mu_u_0*N_i + int_F_edge(x3) - mu_r_L*n(x3) + mu_u_L*N_j = 0
%   rather than via the full flux-at-L formula (see below).
%
% SPATIAL DOMAIN (three compartments, left to right)
%   [0,       L1]              – somatodendritic (SD) region of the SOURCE neuron
%   [L1,      L1+L_ais]        – axon initial segment (AIS)  (reduced diffusivity)
%   [L1+L_ais, L1+L_int-L_syn] – axon proper (diffusion + active transport)
%
% SHOOTING METHOD SUMMARY
%   The unknown is B = n(0), the intracellular tau density at the somatic tip.
%   1. Integrate the SD ODE from x=0 to x=L1  (linear diffusion)
%   2. Continue through the AIS ODE with reduced D (still linear)
%   3. Continue through the axon ODE (nonlinear: includes active transport)
%   4. Evaluate f_init(B): the mass-balance residual at x=L
%   5. fsolve drives f_init(B) → 0 to find the correct B
%
% OUTPUTS
%   J_0  – net tau flux at x=0 (soma end).   Positive = net departure from source node
%   J_L  – net tau flux at x=L (synapse end). Positive = net arrival at target node

close all
format long

% =========================================================================
% SECTION 1: Default parameter values
% =========================================================================
gamma2_0   = 0;       % secondary aggregation rate (unused by default)
frac       = 0.92;    % fraction of tau in the diffusing pool (Konsack 2007)
lambda1_0  = 1e-02;   % AIS diffusivity modifier (D_ais = frac*D*lambda1)
lambda2_0  = lambda1_0;
delta      = 50;      % motor-crowding factor for anterograde transport
epsilon    = 25;      % aggregation-velocity coupling factor

% Compartment lengths in micrometers
L1_    = 200;   % somatodendritic length
L2_    = 200;   % postsynaptic length (legacy; not meshed in current model)
L_int_ = 1000;  % total axon + synaptic region length
L_ais_ = 40;    % AIS length
L_syn_ = 40;    % synaptic terminal / cleft length

% ODE solver tolerances
reltol_    = 1e-6;
abstol_    = 1e-6;
fsolvetol_ = 1e-20;

% Physical scaling
len_scale_  = 1e-03; % um -> mm
time_scale_ = 6 * 30 * 24 * (60)^2; % one month in seconds (the natural NTM time unit)
va = 0.7 * len_scale_ * time_scale_; % anterograde transport speed (mm/month)
vr = 0.7 * len_scale_ * time_scale_; % retrograde transport speed
D  = 12  * len_scale_^2 * time_scale_; % diffusivity

beta     = 1e-06 * time_scale_;   % aggregation saturation constant
gamma1_0 = 0.001 * time_scale_;   % aggregation rate gamma1
gamma1_dt = 0e-05 * time_scale_;  % d(gamma1)/dt (0 = constant)
lambda1_dt = 0e-05;               % d(lambda1)/dt (0 = constant)

% Compute the initial condition for tau concentration.
% The fsolve here finds N_0 such that N_0 + M(N_0) = 2e-2 (init_rescale),
% where M(N) = gamma1*N^2/(beta - gamma2*N) is the bound tau at equilibrium.
% This ensures the total tau (soluble + aggregated) starts at 2e-2.
taufun = @(x) 2e-2 - (x + (gamma1_0 * x.^2) ./ (beta - gamma2_0 * x));
options_taufun = optimset('Display', 'off');
init_rescale_n = fsolve(taufun, 0, options_taufun);

% Initial extracellular tau at the two nodes
N1_0_ = init_rescale_n + 1e-7; % source node (seeded)
N2_0_ = 4e-4;                   % target node (low initial tau)

% Boundary flux parameters
F_edge_0  = 1e-08 * time_scale_; % baseline distributed source along axon
mu_r_0    = 1e-7  * time_scale_; % release rate at x=0
mu_r_L    = mu_r_0;              % release rate at x=L
mu_u_0    = 1e-7  * time_scale_; % uptake rate at x=0
mu_u_L    = mu_u_0;              % uptake rate at x=L
F_edge_dt = 0e-05; % d(F_edge)/dt
F_edge_dx = 0e-05; % d(F_edge)/dx

resmesh_ = 'fine';

% =========================================================================
% SECTION 2: inputParser
% =========================================================================
ip = inputParser;
validScalar = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
addParameter(ip, 'beta',      beta,      validScalar);
addParameter(ip, 'gamma1_0',  gamma1_0,  validScalar);
addParameter(ip, 'gamma2_0',  gamma2_0,  validScalar);
addParameter(ip, 'delta',     delta,     validScalar);
addParameter(ip, 'D',         D,         validScalar);
addParameter(ip, 'epsilon',   epsilon,   validScalar);
addParameter(ip, 'F_edge_0',  F_edge_0);
addParameter(ip, 'frac',      frac,      validScalar);
addParameter(ip, 'lambda1_0', lambda1_0, validScalar);
addParameter(ip, 'lambda2_0', lambda2_0, validScalar);
addParameter(ip, 'L_int',     L_int_,    validScalar);
addParameter(ip, 'L1',        L1_,       validScalar);
addParameter(ip, 'L2',        L2_,       validScalar);
addParameter(ip, 'N1_0',      N1_0_,     validScalar);
addParameter(ip, 'N2_0',      N2_0_,     validScalar);
addParameter(ip, 'mu_r_0',    mu_r_0);
addParameter(ip, 'mu_r_L',    mu_r_L);
addParameter(ip, 'mu_u_0',    mu_u_0);
addParameter(ip, 'mu_u_L',    mu_u_L);
addParameter(ip, 'resmesh',   resmesh_);
addParameter(ip, 'L_ais',     L_ais_);
addParameter(ip, 'L_syn',     L_syn_);
addParameter(ip, 'va',        va);
addParameter(ip, 'vr',        vr);
addParameter(ip, 'fsolvetol', fsolvetol_, validScalar);
addParameter(ip, 'reltol',    reltol_,    validScalar);
addParameter(ip, 'abstol',    abstol_,    validScalar);
addParameter(ip, 'len_scale', len_scale_, validScalar);
addParameter(ip, 'time_scale', time_scale_, validScalar);
parse(ip, varargin{:});

% =========================================================================
% SECTION 3: Convert lengths to model units
% =========================================================================
L1_new    = ip.Results.L1    * ip.Results.len_scale;
L_int_new = ip.Results.L_int * ip.Results.len_scale;
L_ais_new = ip.Results.L_ais * ip.Results.len_scale;
L_syn_new = ip.Results.L_syn * ip.Results.len_scale;
L_total   = L1_new + L_int_new - L_syn_new;

% =========================================================================
% SECTION 4: Build the spatial mesh
% Same logic as Flux_Calculator_RU but without the postsynaptic xmesh2 region.
% =========================================================================
if strcmp(ip.Results.resmesh, 'fine')
    num_comp = 1000;
    num_ext  = 100;
    num_int  = num_comp - 2 * (num_ext);
    xmesh1 = [linspace(0, L1_new - 10*ip.Results.len_scale, num_ext-40), ...
               (L1_new - 9.75*ip.Results.len_scale) : 0.25*ip.Results.len_scale : L1_new];
    xmesh_int = [(L1_new + 0.25*ip.Results.len_scale) : 0.25*ip.Results.len_scale : (L1_new + L_ais_new), ...
                  linspace(L1_new + L_ais_new + 0.25*ip.Results.len_scale, L1_new + L_int_new - L_syn_new, ...
                           num_int - ((L_ais_new + L_syn_new) / (0.25*ip.Results.len_scale))), ...
                 (L1_new + L_int_new - (L_syn_new - 0.25*ip.Results.len_scale)) : 0.25*ip.Results.len_scale : (L1_new + L_int_new - L_syn_new)];
elseif strcmp(ip.Results.resmesh, 'coarse')
    num_comp = 250;
    num_ext  = 25;
    num_int  = num_comp - 2 * (num_ext);
    xmesh1 = [linspace(0, L1_new - 10*ip.Results.len_scale, num_ext-5), ...
               (L1_new - 8*ip.Results.len_scale) : 2*ip.Results.len_scale : L1_new];
    xmesh_int = [(L1_new + 2*ip.Results.len_scale) : 2*ip.Results.len_scale : (L1_new + L_ais_new), ...
                  linspace(L1_new + L_ais_new + 2*ip.Results.len_scale, L1_new + L_int_new - L_syn_new, ...
                           num_int - ((L_ais_new + L_syn_new) / (2*ip.Results.len_scale))), ...
                 (L1_new + L_int_new - (L_syn_new - 2*ip.Results.len_scale)) : 2*ip.Results.len_scale : (L1_new + L_int_new - L_syn_new)];
end
xmesh = [xmesh1, xmesh_int];
xmesh(end)   % print final mesh point (for debugging)

% Extract sub-meshes and boundaries
presyn_mask2 = spatial_mask('presyn');
xmesh_presyn = xmesh(presyn_mask2);
x1 = xmesh_presyn(end);    % presyn / AIS boundary

ais_mask2 = spatial_mask('ais');
xmesh_ais = xmesh(ais_mask2);
x2 = xmesh_ais(end);       % AIS / axon boundary

axon_mask2 = spatial_mask('axon');
xmesh_axon = xmesh(axon_mask2);
x3 = xmesh_axon(end);      % synaptic end of axon
x3   % print (debug)

xmesh = [xmesh_presyn xmesh_ais xmesh_axon]; % trim to active domain

% =========================================================================
% SECTION 5: Time-dependent rate functions (quasi-static: t is frozen)
% =========================================================================
gamma1  = @(t)  ip.Results.gamma1_0 + t .* gamma1_dt;
lambda1 = @(t)  ip.Results.lambda1_0 + t .* lambda1_dt;
F_edge  = @(x, t) (F_edge_0 + t .* F_edge_dt + x .* F_edge_dx);

% =========================================================================
% SECTION 6: Solve the shooting problem at t=0 and report fluxes
% =========================================================================
f0 = [0.1];  % initial guess for shooting parameter B = n(0)
t  = 0;      % time_qs (quasi-static time)

tic
[~, res, flux_0, flux_L] = flux_calculator_RU(xmesh, va, vr, D, t, f0);
time = toc;

J_0 = flux_0;
J_L = flux_L;

fprintf('Residual error %e \n', abs(res))
fprintf('Execution time %e \n', time)


% =========================================================================
% NESTED FUNCTION: flux_calculator_RU
%   Performs the shooting method using the Neumann BC formulation.
%
%   The key difference from Flux_Calculator_RU is the shooting residual:
%     f_init(B) = -mu_r_0*B + mu_u_0*N_i + int_F_edge(x3)
%                 - mu_r_L*n(x3) + mu_u_L*N_j
%
%   This is the GLOBAL mass balance: flux into axon at x=0 + production = flux out at x=L.
%   (Equivalent to demanding that the total tau entering the axon equals what exits.)
% =========================================================================
function [A_1, res, flux_0, flux_L] = flux_calculator_RU(xmesh, va, vr, D, t, f0)

    B_1  = N1_0_; % source node extracellular tau
    C_1  = N2_0_; % target node extracellular tau
    idx0 = logical([1]); % indicator: is this edge active? (1 = yes)

    % ------------------------------------------------------------------
    % Cumulative integral of F_edge: int_F_edge(x) = F_edge_0 * (x - x0)
    % (Simplified: F_edge is uniform in x, so integral is just F_edge_0 * length)
    % ------------------------------------------------------------------
    int_F_edge_vec = zeros(1, length(xmesh));
    F_edge_vec = F_edge(xmesh, t);
    for j = 2:length(xmesh)
        int_F_edge_vec(j) = trapz(xmesh(1:j), F_edge_vec(1:j));
    end
    % Switch to the analytic closed form for speed (valid when F_edge_0 is constant)
    x0 = xmesh(1);
    int_F_edge = @(idx, x) (ip.Results.F_edge_0 .* (x - x0)) .* idx;

    % ------------------------------------------------------------------
    % Set up ODEs for each spatial compartment
    % (same structure as Flux_Calculator_RU; see that file for full derivation)
    % ------------------------------------------------------------------
    presyn_mask = spatial_mask('presyn');
    xmesh_presyn = xmesh(presyn_mask);
    n0 = @(A) A;
    options = odeset('RelTol', ip.Results.reltol, 'AbsTol', ip.Results.abstol);
    n_ss_presyn = @(A, idx) ode45(@(x, n) ode_ss_n(x, t, n, A, D, idx), [0, L1_new], n0(A), options);
    n_ss_presyn = @(A, idx, x) deval(n_ss_presyn(A, idx), x);

    ais_mask = spatial_mask('ais');
    xmesh_ais = xmesh(ais_mask);
    n_ss_ais = @(A, idx) ode45(@(x, n) ode_ss_n(x, t, n, A, D*lambda1(t), idx), ...
                                [L1_new, L1_new+L_ais_new], n_ss_presyn(A, idx, x1), options);
    n_ss_ais = @(A, idx, x) deval(n_ss_ais(A, idx), x);

    axon_mask = spatial_mask('axon');
    xmesh_axon = xmesh(axon_mask);
    n_ss_axon = @(A, idx) ode45(@(x, n) ode_ss_axon(x, t, n, A, idx), ...
                                 [L1_new+L_ais_new, L1_new+L_int_new-L_syn_new], ...
                                 n_ss_ais(A, idx, x2), options);
    n_ss_axon = @(A, idx, x) deval(n_ss_axon(A, idx), x);

    % ------------------------------------------------------------------
    % Shooting residual (Neumann / mass-balance form)
    % ------------------------------------------------------------------
    f_init = @(A) -mu_r_0 .* A + mu_u_0 .* B_1 + int_F_edge(idx0, x3) ...
                  - mu_r_L .* n_ss_axon(A, idx0, x3) + mu_u_L .* C_1;

    options = optimset('TolFun', ip.Results.fsolvetol, 'Display', 'off');
    A_1 = fsolve(f_init, f0, options);  % find B = n(0) such that BC is satisfied
    res = f_init(A_1);                  % residual (quality check)

    % Evaluate full profile for plotting
    n_ss_presyn_plt = n_ss_presyn(A_1, idx0, xmesh_presyn);
    n_ss_ais_plt    = n_ss_ais(A_1, idx0, xmesh_ais);
    n_ss_axon_plt   = n_ss_axon(A_1, idx0, xmesh_axon);
    n_plot = [n_ss_presyn_plt n_ss_ais_plt n_ss_axon_plt];
    figure
    plot(xmesh, n_plot)  % visualize n(x) profile along the axon

    % Output fluxes (sign convention: positive = leaving source / arriving at target)
    flux_0 = -mu_r_0 .* A_1 + mu_u_0 .* B_1
    flux_L =  mu_r_L .* n_ss_axon(A_1, idx0, x3) - mu_u_L .* C_1

    % -------  Nested ODEs  -------

    function nprime = ode_ss_n(x, ~, ~, A, D, idx)
        % Presyn / AIS: purely diffusive steady state
        % dn/dx = (1/D)*(mu_r_0*n(0) - mu_u_0*N_i - int_F_edge)
        nprime = 1./D .* (mu_r_0 .* A - mu_u_0 .* B_1 - int_F_edge(idx, x));
    end

    function nprime = ode_ss_axon(x, t, n, A, idx)
        % Axon proper: diffusion + active transport + nonlinear aggregation coupling
        % dn/dx = (1/(frac*D)) * [(1-frac)*v_eff(n)*n + source - uptake]
        nprime = 1 ./ (ip.Results.frac * ip.Results.D) .* ...
                 ((1 - ip.Results.frac) .* ...
                  ((va * (1 + ip.Results.delta .* n) .* ...
                    (1 - (gamma1(t) * ip.Results.epsilon .* n.^2) ./ ...
                         (ip.Results.beta - ip.Results.gamma2_0 .* n)) - vr)) .* n + ...
                  mu_r_0 .* A - mu_u_0 .* B_1 - int_F_edge(idx, x));
    end

end  % flux_calculator_RU


% =========================================================================
% NESTED FUNCTION: spatial_mask
% =========================================================================
function [maskvals] = spatial_mask(compartment)
    switch compartment
        case 'presyn'
            maskvals = (xmesh <= L1_new);
        case 'ais'
            maskvals = logical(-1 + (xmesh > L1_new) + (xmesh < (L1_new + L_ais_new)));
        case 'axon'
            maskvals = logical(-1 + (xmesh >= L1_new + L_ais_new) + ...
                               (xmesh <= (L1_new + L_int_new - L_syn_new)));
        case 'syncleft'
            maskvals = logical(-1 + (xmesh > L1_new + L_int_new - L_syn_new) + ...
                               (xmesh < (L1_new + L_int_new)));
        case 'postsyn'
            maskvals = logical(-1 + (xmesh >= L1_new + L_int_new) + ...
                               (xmesh <= (L1_new + L_int_new + L2_new)));
        otherwise
            error('Incorrect compartment specification')
    end
end

end  % Flux_Calculator_Neu_edge
