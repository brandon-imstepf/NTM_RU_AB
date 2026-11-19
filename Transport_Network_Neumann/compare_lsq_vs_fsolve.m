% Compare two NTM_ru_ab runs made with the same parameters but different
% shooting solvers (lsqnonlin w/ B>=0 bound, vs the original unconstrained
% fsolve) -- sanity-checks negativity and quantifies how much the fix moved
% the actual tau burden.
% Read-only: loads existing .mat outputs, writes figures to out_dir only.
% Does not touch anything inside the NTM_ru_ab repo's model code.

cd('C:\Users\USER\Documents\Code\NTM_ru_ab\Transport_Network_Neumann');
out_dir = 'C:\Users\USER\Documents\Code\NTM_ru_ab\Transport_Network_Neumann\figures';
if ~isfolder(out_dir), mkdir(out_dir); end

% -------------------------------------------------------------------------
% EDIT THESE TWO LINES: point at your lsqnonlin run and your fsolve
% baseline run (same coupling_term/kappa config, different solver).
% -------------------------------------------------------------------------
lsq_file    = 'call_NTM_brandontest_lsq_mu_u_20261105_194021.mat';
fsolve_file = 'call_NTM_brandontest_fsolve_REPLACE_ME.mat';  % <- fill in after rerun

if ~isfile(fsolve_file)
    error(['fsolve baseline file not found: %s\n' ...
           'Run call_NTM_brandontest.m with the fsolve block active first ' ...
           '(see chat for how), then point fsolve_file at the saved .mat.'], fsolve_file);
end

lsq_d    = load(lsq_file);
fsolve_d = load(fsolve_file);

%% Palette (dataviz skill reference, light mode) -- same constants as compare_sims.m
blue   = double([0x2a 0x78 0xd6])/255;  % fsolve (original)
orange = double([0xeb 0x68 0x34])/255;  % lsqnonlin (bounded)
diverge_neg = double([0xd0 0x3b 0x3b])/255;
diverge_pos = double([0x2a 0x78 0xd6])/255;
ink_primary   = double([0x0b 0x0b 0x0b])/255;
ink_secondary = double([0x52 0x51 0x4e])/255;
ink_muted     = double([0x89 0x87 0x81])/255;
grid_hair     = double([0xe1 0xe0 0xd9])/255;
surface_c     = double([0xfc 0xfc 0xfb])/255;

%% Sanity checks: did either run produce unphysical values?
fprintf('--- Negativity / NaN sanity check ---\n');
fprintf('  fsolve:    min(N)=%.6g  min(M)=%.6g  any NaN/Inf: %d\n', ...
    min(fsolve_d.N(:)), min(fsolve_d.M(:)), any(~isfinite([fsolve_d.N(:); fsolve_d.M(:)])));
fprintf('  lsqnonlin: min(N)=%.6g  min(M)=%.6g  any NaN/Inf: %d\n', ...
    min(lsq_d.N(:)), min(lsq_d.M(:)), any(~isfinite([lsq_d.N(:); lsq_d.M(:)])));

% Time axes must match to compare pointwise
t = lsq_d.t;
assert(isequal(size(t), size(fsolve_d.t)), 'Time grids have different lengths -- not directly comparable.');
if ~isequal(t, fsolve_d.t)
    warning('Time values differ slightly between runs; using the lsqnonlin run''s t for plotting.');
end
nregions = size(lsq_d.N, 1);

total_lsq    = lsq_d.N    + lsq_d.M;   % [nregions x nt]
total_fsolve = fsolve_d.N + fsolve_d.M;

%% Figure 1: global tau burden over time, both solvers
fig1 = figure('Color','w','Position',[100 100 900 550]);
ax = axes(); hold(ax,'on');
plot(t, sum(total_fsolve,1), 'Color', blue,   'LineWidth', 2.5, 'DisplayName', 'fsolve (original)');
plot(t, sum(total_lsq,1),    'Color', orange, 'LineWidth', 2.5, 'DisplayName', 'lsqnonlin (B \geq 0)');
hold(ax,'off');
set(ax,'Color',surface_c,'XColor',ink_muted,'YColor',ink_muted,'GridColor',grid_hair);
grid(ax,'on');
xlabel('time'); ylabel('global total tau (summed over regions)');
title('Global tau burden: fsolve vs lsqnonlin, same parameters', 'Color',ink_primary,'FontWeight','normal');
lgd = legend('Location','best'); lgd.Box = 'off'; lgd.TextColor = ink_secondary;
exportgraphics(fig1, fullfile(out_dir,'solver_compare_fig1_global_burden.png'), 'Resolution',150);
fprintf('Saved solver_compare_fig1_global_burden.png\n');

%% Figure 2: per-region relative difference at final time
fig2 = figure('Color','w','Position',[100 100 900 500]);
final_lsq    = total_lsq(:,end);
final_fsolve = total_fsolve(:,end);
rel_diff = (final_lsq - final_fsolve) ./ max(abs(final_fsolve), eps);
[sorted_diff, order] = sort(rel_diff, 'descend');
bar_colors = zeros(nregions,3);
bar_colors(sorted_diff>=0,:) = repmat(diverge_pos, sum(sorted_diff>=0), 1);
bar_colors(sorted_diff<0,:)  = repmat(diverge_neg, sum(sorted_diff<0), 1);
b = bar(sorted_diff, 'FaceColor','flat');
b.CData = bar_colors;
set(gca,'Color',surface_c,'XColor',ink_muted,'YColor',ink_muted);
xlabel('region (sorted by relative effect size)');
ylabel('relative \Delta total tau at final t  (lsqnonlin - fsolve) / fsolve');
title('Per-region deviation introduced by the B \geq 0 bound', 'Color',ink_primary,'FontWeight','normal');
xticks(1:nregions); xticklabels(string(order)); xtickangle(90);
exportgraphics(fig2, fullfile(out_dir,'solver_compare_fig2_per_region_reldiff.png'), 'Resolution',150);
fprintf('Saved solver_compare_fig2_per_region_reldiff.png\n');

fprintf('\nMax |relative diff| at final time: %.4g (region %d)\n', ...
    max(abs(rel_diff)), order(find(abs(sorted_diff) == max(abs(rel_diff)), 1)));
fprintf('All figures written to %s\n', out_dir);
