%time_qs=0.01;
function [J_0, J_L] = Flux_Calculator_RU(time_qs, varargin)
% Flux_Calculator_RU  Standalone diagnostic: solve the release-uptake edge PDE.
%
% PURPOSE
%   Computes the steady-state flux of intracellular tau at both ends of a
%   SINGLE axon (edge) connecting two brain-region nodes.  This is an
%   earlier, self-contained version of the edge solver; the production-ready
%   version embedded in NetworkFluxCalculator_Neu replaces it.
%
% BIOLOGICAL CONTEXT
%   Each axon is modeled as a 1-D domain [0, L_total] with three spatial
%   compartments:
%     [0,       L1]              – presynaptic somatodendritic (SD) region
%                                  (soma + dendrites of the SOURCE neuron)
%     [L1,      L1+L_ais]        – axon initial segment (AIS): same ODE as
%                                  the SD region but diffusivity is scaled
%                                  by lambda1 (AIS acts as a partial barrier)
%     [L1+L_ais, L1+L_int-L_syn] – axon proper: full active transport ODE
%
%   The variable n(x) is the density of INTRACELLULAR mobile tau along the axon.
%   At steady state, n satisfies different ODEs in each compartment (see below).
%
% BOUNDARY CONDITIONS (Neumann / flux-type)
%   At x=0 (soma):
%     Release INTO N_i (extracellular at source):  mu_r_0 * n(0)
%     Uptake FROM N_i into axon:                   mu_u_0 * N1_0
%     Steady-state BC: D*dn/dx|_0 = mu_r_0*n(0) - mu_u_0*N1_0
%
%   At x=L (synapse):
%     Release FROM axon INTO N_j (extracellular at target): mu_r_L * n(L)
%     Uptake FROM N_j into axon:                             mu_u_L * N2_0
%     Steady-state BC: flux = mu_r_L*n(L) - mu_u_L*N2_0
%
% SHOOTING METHOD
%   n(0) = B is treated as an unknown (the "shooting parameter").
%   The function integrates the steady-state ODE from x=0 to x=L,
%   then solves for B via fsolve so that the global mass-balance condition
%   (f_init = 0, see below) is satisfied.
%
% OUTPUTS
%   J_0 – net flux at x=0: = -mu_r_0*n(0) + mu_u_0*N1_0
%         Positive J_0  → net tau flow FROM the source node INTO the axon
%   J_L – net flux at x=L: = mu_r_L*n(L) - mu_u_L*N2_0
%         Positive J_L  → net tau flow FROM the axon INTO the target node
%
% INPUT
%   time_qs – quasi-static time (used only if gamma1 or lambda1 have
%             time-dependent terms; in the default (dt=0) case, ignored)
%   varargin – optional name-value pairs to override defaults (see inputParser)

close all
format long

% =========================================================================
% SECTION 1: Default parameter values
% All lengths are in micrometers (um), velocities in um/s, times in seconds
% unless explicitly rescaled by len_scale and time_scale below.
% =========================================================================
gamma2_0   = 0;         % secondary aggregation rate (usually set to 0)
frac       = 0.92;      % fraction of tau in diffusing (free) pool (Konsack 2007)
lambda1_0  = 1e-02;     % AIS diffusivity factor (D_ais = D * lambda1)
lambda2_0  = lambda1_0; % (legacy; matches lambda1 here)
delta      = 1;         % motor-crowding enhancement (v_a scales as 1 + delta*n)
epsilon    = 1e-02;     % aggregation-velocity coupling

% Compartment lengths (um, before scaling)
L1_   = 200;    % somatodendritic region length
L2_   = 200;    % postsynaptic region length (legacy; L2 not in current mesh)
L_int_ = 1000;  % total axon + synapse length
L_ais_ = 40;    % axon initial segment length
L_syn_ = 40;    % synaptic cleft / terminal length

% Initial extracellular tau concentrations at the two nodes
N1_0_ = 1e-03;  % extracellular tau at source node (x=0 boundary)
N2_0_ = 5e-03;  % extracellular tau at target node (x=L boundary)
M1_0_ = 0;      % bound tau at source (not used in PDE but kept for mass balance)
M2_0_ = 0;      % bound tau at target (not used in PDE but kept for mass balance)

% ODE solver tolerances
reltol_    = 1e-10;
abstol_    = 1e-10;
fsolvetol_ = 1e-20;

% Physical scaling factors
% In the NTM the natural time unit is one month (~2.6e6 s) and length unit is mm.
len_scale_  = 1e-03;  % convert um -> mm  (all lengths multiplied by this)
time_scale_ = 1;      % set to 1 here (no time rescaling in this standalone script)

% Biophysical constants (Konsack 2007, scaled to model units)
va = 0.7  * len_scale_ * time_scale_;   % anterograde active transport velocity
vr = 0.7  * len_scale_ * time_scale_;   % retrograde active transport velocity
D  = 12   * len_scale_^2 * time_scale_; % diffusivity of free tau

% Aggregation / kinetic constants
beta     = 1e-06 * time_scale_; % aggregation saturation constant
gamma1_0 = 1e-05 * time_scale_; % aggregation rate (soluble -> bound)
gamma1_dt = 0e-05 * time_scale_; % time derivative of gamma1 (0 = time-invariant)
lambda1_dt = 0e-05;              % time derivative of lambda1 (0 = time-invariant)

% Boundary flux parameters
F_edge_0  = 1e-08 * time_scale_; % baseline distributed tau source along axon
mu_r_0    = 1e-05;               % release rate at soma end (x=0)
mu_r_L    = mu_r_0;              % release rate at synapse end (x=L)
mu_u_0    = 1e-04;               % uptake rate at soma end
mu_u_L    = mu_u_0;              % uptake rate at synapse end
F_edge_dt = 0e-05;               % time derivative of F_edge (0 = constant)
F_edge_dx = 0e-05;               % spatial gradient of F_edge (0 = uniform)

resmesh_    = 'fine';
total_mass_ = 184;   % legacy total mass target (not actively used)

% =========================================================================
% SECTION 2: inputParser – override any default with name-value pair args
% =========================================================================
ip = inputParser;
validScalar = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
addParameter(ip, 'beta',       beta,      validScalar);
addParameter(ip, 'gamma1_0',   gamma1_0,  validScalar);
addParameter(ip, 'gamma2_0',   gamma2_0,  validScalar);
addParameter(ip, 'delta',      delta,     validScalar);
addParameter(ip, 'D',          D,         validScalar);
addParameter(ip, 'epsilon',    epsilon,   validScalar);
addParameter(ip, 'F_edge_0',   F_edge_0);
addParameter(ip, 'frac',       frac,      validScalar);
addParameter(ip, 'lambda1_0',  lambda1_0, validScalar);
addParameter(ip, 'lambda2_0',  lambda2_0, validScalar);
addParameter(ip, 'L_int',      L_int_,    validScalar);
addParameter(ip, 'L1',         L1_,       validScalar);
addParameter(ip, 'L2',         L2_,       validScalar);
addParameter(ip, 'N1_0',       N1_0_,     validScalar);
addParameter(ip, 'N2_0',       N2_0_,     validScalar);
addParameter(ip, 'M1_0',       M1_0_,     validScalar);
addParameter(ip, 'M2_0',       M2_0_,     validScalar);
addParameter(ip, 'mu_r_0',     mu_r_0);
addParameter(ip, 'mu_r_L',     mu_r_L);
addParameter(ip, 'mu_u_0',     mu_u_0);
addParameter(ip, 'mu_u_L',     mu_u_L);
addParameter(ip, 'resmesh',    resmesh_);
addParameter(ip, 'L_ais',      L_ais_);
addParameter(ip, 'L_syn',      L_syn_);
addParameter(ip, 'total_mass', total_mass_);
addParameter(ip, 'va',         va);
addParameter(ip, 'vr',         vr);
addParameter(ip, 'fsolvetol',  fsolvetol_, validScalar);
addParameter(ip, 'reltol',     reltol_,    validScalar);
addParameter(ip, 'abstol',     abstol_,    validScalar);
addParameter(ip, 'len_scale',  len_scale_, validScalar);
addParameter(ip, 'time_scale', time_scale_, validScalar);
parse(ip, varargin{:});

% =========================================================================
% SECTION 3: Convert lengths from um to model units (mm by default)
% =========================================================================
L1_new    = ip.Results.L1    * ip.Results.len_scale;
L2_new    = ip.Results.L2    * ip.Results.len_scale;
L_int_new = ip.Results.L_int * ip.Results.len_scale;
L_ais_new = ip.Results.L_ais * ip.Results.len_scale;
L_syn_new = ip.Results.L_syn * ip.Results.len_scale;

% Total spatial domain: soma + axon (minus the synaptic terminal)
L_total = L1_new + L2_new + L_int_new;

% =========================================================================
% SECTION 4: Build the spatial mesh (xmesh)
%
% Two resolution options:
%   'fine'   – 1000 points total, denser near the compartment boundaries
%   'coarse' – 250 points total, suitable for fast parameter sweeps
%
% The mesh covers three regions:
%   xmesh1    – somatodendritic compartment [0, L1]
%   xmesh_int – AIS + axon proper [L1, L1+L_int-L_syn]
%   xmesh2    – postsynaptic region [L1+L_int, L_total]  (legacy; not used in current model)
% =========================================================================
if strcmp(ip.Results.resmesh, 'fine')
    num_comp = 1000;   % total mesh points
    num_ext  = 200;    % points per external (SD) compartment
    num_int  = num_comp - 2 * (num_ext);
    % Fine uniform grid for the SD region, with extra density near boundary L1
    xmesh1 = [linspace(0, L1_new - 10*ip.Results.len_scale, num_ext-40), ...
               (L1_new - 9.75*ip.Results.len_scale) : 0.25*ip.Results.len_scale : L1_new];
    % Postsynaptic mesh (symmetric to soma)
    xmesh2 = [(L1_new + L_int_new + 0.25*ip.Results.len_scale) : 0.25*ip.Results.len_scale : ...
               (L1_new + L_int_new + 10*ip.Results.len_scale), ...
               linspace(L1_new + L_int_new + 10.25*ip.Results.len_scale, L_total, num_ext-40)];
    % Interior (AIS + axon) mesh with fine spacing at AIS/synapse transitions
    xmesh_int = [(L1_new + 0.25*ip.Results.len_scale) : 0.25*ip.Results.len_scale : (L1_new + L_ais_new), ...
                  linspace(L1_new + L_ais_new + 0.25*ip.Results.len_scale, L1_new + L_int_new - L_syn_new, ...
                           num_int - ((L_ais_new + L_syn_new) / (0.25*ip.Results.len_scale))), ...
                 (L1_new + L_int_new - (L_syn_new - 0.25*ip.Results.len_scale)) : 0.25*ip.Results.len_scale : (L1_new + L_int_new)];

elseif strcmp(ip.Results.resmesh, 'coarse')
    num_comp = 250;
    num_ext  = 25;
    num_int  = num_comp - 2 * (num_ext);
    xmesh1 = [linspace(0, L1_new - 10*ip.Results.len_scale, num_ext-5), ...
               (L1_new - 8*ip.Results.len_scale) : 2*ip.Results.len_scale : L1_new];
    xmesh2 = [(L1_new + L_int_new + 2*ip.Results.len_scale) : 2*ip.Results.len_scale : ...
               (L1_new + L_int_new + 10*ip.Results.len_scale), ...
               linspace(L1_new + L_int_new + 12*ip.Results.len_scale, L_total, num_ext-5)];
    xmesh_int = [(L1_new + 2*ip.Results.len_scale) : 2*ip.Results.len_scale : (L1_new + L_ais_new), ...
                  linspace(L1_new + L_ais_new + 2*ip.Results.len_scale, L1_new + L_int_new - L_syn_new, ...
                           num_int - ((L_ais_new + L_syn_new) / (2*ip.Results.len_scale))), ...
                 (L1_new + L_int_new - (L_syn_new - 2*ip.Results.len_scale)) : 2*ip.Results.len_scale : (L1_new + L_int_new)];
end
xmesh = [xmesh1, xmesh_int, xmesh2];

% Extract sub-meshes for each compartment
presyn_mask2 = spatial_mask('presyn');
xmesh_presyn = xmesh(presyn_mask2);
x1 = xmesh_presyn(end);   % boundary between presyn and AIS

ais_mask2 = spatial_mask('ais');
xmesh_ais = xmesh(ais_mask2);
x2 = xmesh_ais(end);      % boundary between AIS and axon

axon_mask2 = spatial_mask('axon');
xmesh_axon = xmesh(axon_mask2);
x3 = xmesh_axon(end);     % synaptic end of axon

% Trim the combined mesh to only the active compartments (presyn + AIS + axon)
xmesh = [xmesh_presyn xmesh_ais xmesh_axon];

% =========================================================================
% SECTION 5: Time-dependent coefficient functions
% In the quasi-static steady-state, t is fixed (= time_qs), so these
% effectively return constants unless gamma1_dt or lambda1_dt are nonzero.
% =========================================================================
gamma1 = @(t) ip.Results.gamma1_0 + t .* gamma1_dt;  % aggregation rate at time t
lambda1 = @(t) ip.Results.lambda1_0 + t .* lambda1_dt; % AIS permeability at time t

% Distributed tau production along the axon (constant by default)
F_edge = @(x, t) F_edge_0 + t .* F_edge_dt + x .* F_edge_dx;

% =========================================================================
% SECTION 6: Solve for the shooting parameter and compute fluxes
% =========================================================================
f0 = 0;       % initial guess for the shooting parameter B = n(0)
t  = time_qs; % the quasi-static time point

tic
[~, ~, res, flux_0, flux_L] = flux_calculator_RU(xmesh, va, vr, D, t, f0);
time = toc;

J_0 = flux_0;
J_L = flux_L;

fprintf('Residual error %e \n', abs(res))
fprintf('Execution time %e \n', time)


% =========================================================================
% NESTED FUNCTION: flux_calculator_RU
%   Performs the shooting solve and returns the axon fluxes.
%   n_plot is the full spatial profile of n(x) for visualization.
% =========================================================================
function [n_plot, A_1, res, flux_0, flux_L] = flux_calculator_RU(xmesh, va, vr, D, t, f0)

    B_1 = N1_0_; % extracellular tau concentration at source node (x=0)
    C_1 = N2_0_; % extracellular tau concentration at target node (x=L)

    % ------------------------------------------------------------------
    % Pre-compute the cumulative integral of the distributed source F_edge
    % along xmesh.  int_F_edge(x) = integral from 0 to x of F_edge(x',t) dx'.
    % This enters the ODE as a cumulative production term.
    % ------------------------------------------------------------------
    int_F_edge = zeros(1, length(xmesh));
    F_edge_vec = F_edge(xmesh, t);
    for j = 2:length(xmesh)
        int_F_edge(j) = trapz(xmesh(1:j), F_edge_vec(1:j));
    end
    int_F_edge = @(x) interp1(xmesh, int_F_edge, x); % callable at any x

    % ------------------------------------------------------------------
    % Steady-state ODE in the somatodendritic region [0, L1]
    %
    % In this purely diffusive compartment, the steady-state reads:
    %   frac * D * d²n/dx² = - (source - uptake - production)
    %
    % In 1-D with constant coefficients, this simplifies to:
    %   dn/dx = (1/D) * (mu_r_0 * A - mu_u_0 * B_1 - int_F_edge(x))
    %
    % where A = n(0) is the shooting parameter (initial condition for ODE).
    % Note: n(0) = A is both the IC and appears in the RHS because release
    % (mu_r_0 * n(0)) is treated as a constant source set at the somatic value.
    % ------------------------------------------------------------------
    presyn_mask = spatial_mask('presyn');
    xmesh_presyn = xmesh(presyn_mask);
    n0 = @(A) A;   % IC: n at x=0 equals the shooting parameter
    options = odeset('RelTol', ip.Results.reltol, 'AbsTol', ip.Results.abstol, ...
                     'NonNegative', 1:length(n0));
    n_ss_presyn = @(A) ode15s(@(x, n) ode_ss_n(x, t, n, A, D), [0, L1_new], n0(A), options);
    n_ss_presyn = @(A, x) deval(n_ss_presyn(A), x); % evaluate at arbitrary x

    % ------------------------------------------------------------------
    % Steady-state ODE in the AIS [L1, L1+L_ais]
    %
    % Same structure as the presyn ODE but diffusivity is reduced by lambda1:
    %   D_ais = D * lambda1(t)
    % This models the AIS as a partial diffusion barrier.
    % IC: continuity of n at x = x1 (end of presyn region).
    % ------------------------------------------------------------------
    ais_mask = spatial_mask('ais');
    xmesh_ais = xmesh(ais_mask);
    n0 = @(A) n_ss_presyn(A, x1);  % continuity BC at presyn/AIS boundary
    n_ss_ais = @(A) ode15s(@(x, n) ode_ss_n(x, t, n, A, D*lambda1(t)), ...
                            [L1_new, L1_new+L_ais_new], n0(A), options);
    n_ss_ais = @(A, x) deval(n_ss_ais(A), x);

    % ------------------------------------------------------------------
    % Steady-state ODE in the axon proper [L1+L_ais, L1+L_int-L_syn]
    %
    % In the axon, tau undergoes ACTIVE TRANSPORT (not just diffusion).
    % The effective velocity depends on n through crowding (delta) and
    % aggregation (gamma1, epsilon):
    %
    %   v_eff = v_a * (1 + delta*n) * (1 - gamma1*epsilon*n^2/(beta-gamma2*n)) - v_r
    %
    %   v_a*(1 + delta*n)         – anterograde velocity, enhanced when n is high
    %                               (crowded motors get recruited; Konsack model)
    %   *(1 - gamma1*epsilon*n^2/(beta - gamma2*n))
    %                             – reduction because aggregated tau (M ∝ n^2)
    %                               impedes motor function
    %   - v_r                     – subtracted retrograde velocity
    %
    % The full steady-state ODE integrates diffusion + transport + source/uptake:
    %   dn/dx = (1/(frac*D)) * [(1-frac)*v_eff*n + mu_r_0*A - mu_u_0*B1 - F_edge(x)]
    %
    % IC: continuity at the AIS/axon boundary x = x2.
    % ------------------------------------------------------------------
    axon_mask = spatial_mask('axon');
    xmesh_axon = xmesh(axon_mask);
    n0 = @(A) n_ss_ais(A, x2);   % continuity BC at AIS/axon boundary
    n_ss_axon = @(A) ode15s(@(x, n) ode_ss_axon(x, t, n, A), ...
                             [L1_new+L_ais_new, L1_new+L_int_new-L_syn_new], n0(A), options);
    n_ss_axon = @(A, x) deval(n_ss_axon(A), x);

    % ------------------------------------------------------------------
    % Define the shooting residual f_init(A).
    %
    % The global mass-balance condition across the entire axon is:
    %   (net flux in at x=0) + (cumulative production) = (net flux out at x=L)
    %
    % Written out:
    %   f_init(A) = -D*frac*(dn/dx|_{x3}) + (1-frac)*n(x3)*v_eff(x3) - mu_r_L*n(x3) + mu_u_L*C_1
    %
    % This condition is exactly the natural Neumann BC at x = L (the synaptic end).
    % fsolve finds A such that f_init(A) = 0.
    % ------------------------------------------------------------------
    v = @(A, x) ((va * (1 + ip.Results.delta .* n_ss_axon(A, x)) .* ...
                  (1 - gamma1(t) * ip.Results.epsilon .* n_ss_axon(A, x).^2 ./ ...
                  (ip.Results.beta - ip.Results.gamma2_0 .* n_ss_axon(A, x))) - vr));
    f_init = @(A) -ip.Results.D * ip.Results.frac .* ode_ss_axon(x3, t, n_ss_axon(A, x3), A) + ...
                   (1 - ip.Results.frac) .* n_ss_axon(A, x3) .* v(A, x3) - ...
                   mu_r_L * n_ss_axon(A, x3) + mu_u_L * C_1;

    options = optimset('TolFun', ip.Results.fsolvetol, 'Display', 'off');
    A_1 = fsolve(f_init, f0, options);   % solve for the shooting parameter
    res = f_init(A_1);                   % residual (should be ~0 if converged)

    % Evaluate the full spatial profile for plotting
    n_ss_presyn_plt = n_ss_presyn(A_1, xmesh_presyn);
    n_ss_ais_plt    = n_ss_ais(A_1, xmesh_ais);
    n_ss_axon_plt   = n_ss_axon(A_1, xmesh_axon);
    n_plot = [n_ss_presyn_plt n_ss_ais_plt n_ss_axon_plt];

    % ------------------------------------------------------------------
    % Compute the boundary fluxes:
    %   flux_0 = net flux from source node INTO the axon at x=0
    %            = mu_u_0*N_source - mu_r_0*n(0)
    %            (positive → node loses tau to axon)
    %   flux_L = net flux from axon INTO the target node at x=L
    %            = mu_r_L*n(L) - mu_u_L*N_target
    %            (positive → target node gains tau from axon)
    % ------------------------------------------------------------------
    flux_0 = -mu_r_0 .* A_1 + mu_u_0 .* B_1
    flux_L = mu_r_L .* n_ss_axon(A_1, x3) - mu_u_L .* C_1

    % -------  Nested ODEs  -------

    % ODE for presyn and AIS compartments (pure diffusion with source/uptake)
    function nprime = ode_ss_n(x, ~, ~, A, D)
        % dn/dx = (1/D) * (release - uptake - production)
        nprime = 1./D .* (mu_r_0 .* A - mu_u_0 .* B_1 - int_F_edge(x));
    end

    % ODE for axon proper (diffusion + nonlinear active transport)
    function nprime = ode_ss_axon(x, t, n, A)
        % dn/dx = (1/(frac*D)) * [(1-frac)*v_eff*n + release - uptake - production]
        nprime = 1 ./ (ip.Results.frac * ip.Results.D) .* ...
                 ((1 - ip.Results.frac) .* ...
                  ((va * (1 + ip.Results.delta .* n) .* ...
                    (1 - (gamma1(t) * ip.Results.epsilon .* n.^2) ./ ...
                         (ip.Results.beta - ip.Results.gamma2_0 .* n)) - vr)) .* n + ...
                  mu_r_0 .* A - mu_u_0 .* B_1 - int_F_edge(x));
    end
end  % flux_calculator_RU


% =========================================================================
% NESTED FUNCTION: spatial_mask
%   Returns a logical index into xmesh for the named compartment.
% =========================================================================
function [maskvals] = spatial_mask(compartment)
    switch compartment
        case 'presyn'
            % Somatodendritic: x in [0, L1]
            maskvals = (xmesh <= L1_new);
        case 'ais'
            % Axon initial segment: x in (L1, L1+L_ais)
            maskvals = logical(-1 + (xmesh > L1_new) + (xmesh < (L1_new + L_ais_new)));
        case 'axon'
            % Axon proper: x in [L1+L_ais, L1+L_int-L_syn]
            maskvals = logical(-1 + (xmesh >= L1_new + L_ais_new) + ...
                               (xmesh <= (L1_new + L_int_new - L_syn_new)));
        case 'syncleft'
            % Synaptic cleft: x in (L1+L_int-L_syn, L1+L_int)  [not used in current model]
            maskvals = logical(-1 + (xmesh > L1_new + L_int_new - L_syn_new) + ...
                               (xmesh < (L1_new + L_int_new)));
        case 'postsyn'
            % Postsynaptic: x in [L1+L_int, L1+L_int+L2]  [not used in current model]
            maskvals = logical(-1 + (xmesh >= L1_new + L_int_new) + ...
                               (xmesh <= (L1_new + L_int_new + L2_new)));
        otherwise
            error('Incorrect compartment specification')
    end
end

end  % Flux_Calculator_RU
