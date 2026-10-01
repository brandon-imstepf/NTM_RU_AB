function plot_sweep_failures(mat_file)
% PLOT_SWEEP_FAILURES  Visualize which sweep combinations failed to converge.
%
% Loads a fixed_sweep_output_task##.mat file (spread_final/spread_t50/
% burden_final are NaN exactly where a run failed -- see call_NTM_brandontest_sweep.m)
% and draws one heatmap per coupling term: kappa_frac on the y-axis, the 16
% (r_frac, Dif_frac) combinations flattened along the x-axis, red = failed.
%
% Usage: plot_sweep_failures('fixed_sweep_output_task00.mat')

if nargin < 1, mat_file = 'fixed_sweep_output_task00.mat'; end
d = load(mat_file);

fail = isnan(d.spread_final);  % logical, same shape as sz: [term, kappa, r, Dif]
nTerms = numel(d.coupling_terms);
nKappa = numel(d.kappa_frac);
nR     = numel(d.r_frac);
nDif   = numel(d.Dif_frac);

figure('Color', 'w', 'Position', [100 100 1100 300]);
t = tiledlayout(1, nTerms, 'TileSpacing', 'compact', 'Padding', 'compact');

% simple 2-color map: pass = light neutral, fail = red
cmap = [0.95 0.95 0.93; 0.82 0.23 0.23];

for ti = 1:nTerms
    nexttile;
    % collapse (r_frac, Dif_frac) into one axis so each term is a flat 2-D image
    img = reshape(fail(ti,:,:,:), [nKappa, nR*nDif]);
    imagesc(img);
    colormap(gca, cmap);
    clim([0 1]);

    yticks(1:nKappa);
    yticklabels(string(d.kappa_frac));
    ylabel('\kappa_{frac}');

    % x tick per (r_frac, Dif_frac) combo, shown compactly
    xticks(1:nR*nDif);
    [rG, dG] = ndgrid(d.r_frac, d.Dif_frac);
    xlabels = compose('%.1g,%.1g', rG(:), dG(:));
    xticklabels(xlabels);
    xtickangle(90);
    xlabel('(r_{frac}, Dif_{frac})');

    title(d.coupling_terms{ti}, 'Interpreter', 'none');
    set(gca, 'FontSize', 12, 'TickDir', 'out');
    axis tight;
end

n_fail = sum(fail(:));
n_total = numel(fail);
title(t, sprintf('Sweep failures: %d/%d (%.1f%%)', n_fail, n_total, 100*n_fail/n_total));

end
