% Figures illustrating the edge-3 axon blow-up: fixed n* threshold + runaway.
% Read-only: loads the existing debug dump, writes figures to out_dir only.

cd('C:\Users\USER\Documents\Code\NTM_ru_ab\Transport_Network_Neumann');
out_dir = 'C:\Users\USER\Documents\Code\NTM_ru_ab\figures\axon_crash';
if ~isfolder(out_dir), mkdir(out_dir); end

d = load('axon_debug_log_region3_20261026_104500.mat');
term_log = d.term_log;
crashed_region = d.i;
region_ids = arrayfun(@(s) s.region_i, term_log);
rows = term_log(region_ids == crashed_region);
nedges = numel(rows(1).n);

calls = [rows.call];
[calls_sorted, order] = sort(calls);
N_mat = cell2mat(arrayfun(@(s) s.n(:)', rows, 'UniformOutput', false)');
R_mat = cell2mat(arrayfun(@(s) s.R(:)', rows, 'UniformOutput', false)');
N_mat = N_mat(order, :);
R_mat = R_mat(order, :);

% Critical n: R = (gamma1/beta)*epsilon*n^2 = 1  =>  n* = sqrt(1/coef)
% Back out the coefficient directly from the logged data (R and n are both
% logged, so this needs no assumed constants): coef = R ./ n.^2 (use a safe,
% well-away-from-zero sample so the ratio is numerically clean).
safe_row = 1; safe_edge = 1;
coef = R_mat(safe_row, safe_edge) / N_mat(safe_row, safe_edge)^2;
n_star = sqrt(1/coef);
fprintf('Recovered R = %.6g * n^2 (coef), n* (R=1) = %.6g\n', coef, n_star);

%% Legibility defaults -- bigger type throughout
set(groot, 'DefaultAxesFontSize', 13);
set(groot, 'DefaultTextFontSize', 13);
set(groot, 'DefaultAxesLineWidth', 1);

%% Palette (dataviz skill reference, light mode) -- categorical slots 1-6
col = { double([0x2a 0x78 0xd6])/255, ...  % 1 blue
        double([0xeb 0x68 0x34])/255, ...  % 2 orange
        double([0x1b 0xaf 0x7a])/255, ...  % 3 aqua
        double([0xed 0xa1 0x00])/255, ...  % 4 yellow
        double([0xe8 0x7b 0xa4])/255, ...  % 5 magenta
        double([0x62 0x50 0xd6])/255 };    % 6 (reusing violet; 6-edge categorical set)
ink_primary = double([0x0b 0x0b 0x0b])/255;
ink_muted   = double([0x89 0x87 0x81])/255;
surface_c   = double([0xfc 0xfc 0xfb])/255;
grid_hair   = double([0xe1 0xe0 0xd9])/255;

%% Figure 1: |n| (log scale) for all 6 edges vs RHS call index, threshold line
fig1 = figure('Color','w','Position',[100 100 1200 650]);
ax = axes(); hold(ax,'on');
for e = 1:nedges
    lw = 2.5; if e==3 || e==4, lw = 3; end  % the two edges that matter stand out slightly
    semilogy(calls_sorted, abs(N_mat(:,e)), 'Color', col{e}, 'LineWidth', lw, ...
        'DisplayName', sprintf('edge %d', e));
end
yline(n_star, '--', 'Color', ink_muted, 'LineWidth', 1.5, 'HandleVisibility','off');
hold(ax,'off');
set(ax,'YScale','log','Color',surface_c,'XColor',ink_muted,'YColor',ink_muted,'GridColor',grid_hair,'FontSize',13);
grid(ax,'on');
text(ax, 0.015, 0.95, sprintf('dashed line: n^* = %.4g (crowding saturates)', n_star), ...
    'Units','normalized', 'Color',ink_primary, 'FontSize',12, 'BackgroundColor',[1 1 1 0.75], 'Margin',3);
xlabel('RHS call index (solver evaluation order)','FontSize',14); ylabel('|n|  (log scale)','FontSize',14);
title('Only edges 3 and 4 cross the threshold -- edge 3 diverges', ...
    'Color',ink_primary,'FontWeight','bold','FontSize',16);
lgd = legend('Location','southoutside','Orientation','horizontal','NumColumns',6,'FontSize',12);
lgd.Box='off'; lgd.TextColor=ink_primary;
exportgraphics(fig1, fullfile(out_dir,'fig1_all_edges_n_vs_threshold.png'), 'Resolution',150);
fprintf('Saved fig1_all_edges_n_vs_threshold.png\n');

%% Figure 2: R vs n (log-log), theoretical curve + each edge's observed range
fig2 = figure('Color','w','Position',[100 100 950 700]);
ax = axes(); hold(ax,'on');
n_theory = logspace(-4, log10(max(abs(N_mat(:)))*1.2), 400);
loglog(n_theory, coef*n_theory.^2, '-', 'Color', grid_hair, 'LineWidth', 4, 'HandleVisibility','off');
for e = 1:nedges
    ms = 14; if e==3 || e==4, ms = 22; end  % the two edges that matter stand out
    loglog(abs(N_mat(:,e)), R_mat(:,e), '.', 'Color', col{e}, 'MarkerSize', ms, ...
        'DisplayName', sprintf('edge %d', e));
end
yline(1, ':', 'Color', ink_primary, 'LineWidth', 1.2, 'HandleVisibility','off');
xline(n_star, ':', 'Color', ink_primary, 'LineWidth', 1.2, 'HandleVisibility','off');
hold(ax,'off');
set(ax,'XScale','log','YScale','log','Color',surface_c,'XColor',ink_muted,'YColor',ink_muted,'GridColor',grid_hair,'FontSize',13);
grid(ax,'on');
text(ax, n_star*1.3, 1e-3, sprintf('n^* = %.4g', n_star), 'Color',ink_primary, 'FontSize',12, 'Rotation',90);
text(ax, 2e-4, 1.4, 'R = 1 (crowding saturates)', 'Color',ink_primary, 'FontSize',12);
xlabel('|n|  (log scale)','FontSize',14); ylabel('R  (log scale)','FontSize',14);
title('Edges 3 and 4 push n far past where R saturates', ...
    'Color',ink_primary,'FontWeight','bold','FontSize',16);
lgd = legend('Location','southoutside','Orientation','horizontal','NumColumns',6,'FontSize',12);
lgd.Box='off'; lgd.TextColor=ink_primary;
exportgraphics(fig2, fullfile(out_dir,'fig2_R_vs_n_threshold.png'), 'Resolution',150);
fprintf('Saved fig2_R_vs_n_threshold.png\n');

%% Figure 3: zoomed-in runaway -- edge 3's final ~1500 calls, n and R stacked
fig3 = figure('Color','w','Position',[100 100 1050 700]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
zoom_n = max(1, numel(calls_sorted)-1499):numel(calls_sorted);

nexttile;
plot(calls_sorted(zoom_n), N_mat(zoom_n,3), 'Color', col{3}, 'LineWidth', 2.5);
set(gca,'Color',surface_c,'XColor',ink_muted,'YColor',ink_muted,'GridColor',grid_hair,'FontSize',13); grid on;
ylabel('n, edge 3','FontSize',14);
title('Edge 3''s last 1500 solver calls before the crash', 'Color',ink_primary,'FontWeight','bold','FontSize',16);

nexttile;
semilogy(calls_sorted(zoom_n), R_mat(zoom_n,3), 'Color', col{2}, 'LineWidth', 2.5);
set(gca,'Color',surface_c,'XColor',ink_muted,'YColor',ink_muted,'GridColor',grid_hair,'FontSize',13); grid on;
xlabel('RHS call index','FontSize',14); ylabel('R  (log scale)','FontSize',14);

exportgraphics(fig3, fullfile(out_dir,'fig3_edge3_runaway_zoom.png'), 'Resolution',150);
fprintf('Saved fig3_edge3_runaway_zoom.png\n');

fprintf('\nAll figures written to %s\n', out_dir);
