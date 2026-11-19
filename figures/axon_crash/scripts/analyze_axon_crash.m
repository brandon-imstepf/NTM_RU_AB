% Read-only inspection of the axon debug term_log dumped on crash.
% Does not touch anything inside the NTM_ru_ab repo's model code.

cd('C:\Users\USER\Documents\Code\NTM_ru_ab\Transport_Network_Neumann');
out_dir = 'C:\Users\USER\Documents\Code\NTM_ru_ab\figures\axon_crash';
if ~isfolder(out_dir), mkdir(out_dir); end

d = load('axon_debug_log_region3_20261026_104500.mat');
term_log = d.term_log;
crashed_region = d.i;
fprintf('Crashed on region_i = %d. Total RHS evaluations logged (all regions up to crash): %d\n', ...
    crashed_region, numel(term_log));

region_ids = arrayfun(@(s) s.region_i, term_log);
rows = term_log(region_ids == crashed_region);
fprintf('Entries belonging to the crashing region (%d): %d\n', crashed_region, numel(rows));

nedges = numel(rows(1).n);
fprintf('Edge-batch size for this region (number of source edges j -> i=%d): %d\n', crashed_region, nedges);

xs = [rows.x];
fprintf('x range covered: [%.6g, %.6g]\n', min(xs), max(xs));
fprintf('Number of DISTINCT x values visited: %d (out of %d calls -- repeats = step retries)\n', ...
    numel(unique(xs)), numel(xs));

% Stack per-edge traces into [ncalls x nedges] matrices
N_mat   = cell2mat(arrayfun(@(s) s.n(:)',           rows, 'UniformOutput', false)');
R_mat   = cell2mat(arrayfun(@(s) s.R(:)',           rows, 'UniformOutput', false)');
Veff_mat= cell2mat(arrayfun(@(s) s.v_eff(:)',       rows, 'UniformOutput', false)');
NP_mat  = cell2mat(arrayfun(@(s) s.nprime(:)',      rows, 'UniformOutput', false)');

fprintf('\nPer-edge summary over the whole crashing-region log:\n');
for e = 1:nedges
    fprintf('  edge %2d: max(n)=%.6g  max(R)=%.6g  max(|v_eff|)=%.6g  max(|nprime|)=%.6g\n', ...
        e, max(abs(N_mat(:,e))), max(R_mat(:,e)), max(abs(Veff_mat(:,e))), max(abs(NP_mat(:,e))));
end

% Identify the edge with the largest peak R (closest to / past the R=1 crowding-saturation point)
[~, worst_edge] = max(max(R_mat,[],1));
fprintf('\nWorst edge by peak R: edge %d (peak R = %.6g)\n', worst_edge, max(R_mat(:,worst_edge)));

% Order by call index (call order == solver's actual evaluation order,
% NOT monotonic in x -- repeats/backtracks are failed step attempts)
calls = [rows.call];
[calls_sorted, order] = sort(calls);
xs_o   = xs(order);
n_o    = N_mat(order, worst_edge);
R_o    = R_mat(order, worst_edge);
veff_o = Veff_mat(order, worst_edge);
np_o   = NP_mat(order, worst_edge);

fprintf('\nLast 10 RHS evaluations before the crash (worst edge = %d):\n', worst_edge);
fprintf('  %8s %12s %12s %12s %12s %12s\n', 'call', 'x', 'n', 'R', 'v_eff', 'nprime');
last10 = max(1, numel(calls_sorted)-9):numel(calls_sorted);
for k = last10
    fprintf('  %8d %12.6g %12.6g %12.6g %12.6g %12.6g\n', ...
        calls_sorted(k), xs_o(k), n_o(k), R_o(k), veff_o(k), np_o(k));
end

%% Plot: n, R, v_eff, nprime vs call index, for the worst edge
blue   = double([0x2a 0x78 0xd6])/255;
orange = double([0xeb 0x68 0x34])/255;
aqua   = double([0x1b 0xaf 0x7a])/255;
yellow = double([0xed 0xa1 0x00])/255;
ink_primary   = double([0x0b 0x0b 0x0b])/255;
ink_muted     = double([0x89 0x87 0x81])/255;
surface_c     = double([0xfc 0xfc 0xfb])/255;

fig = figure('Color','w','Position',[100 100 1400 800]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
plot(calls_sorted, n_o, 'Color', blue, 'LineWidth', 1.5);
set(gca,'Color',surface_c,'XColor',ink_muted,'YColor',ink_muted);
xlabel('RHS call index'); ylabel('n');
title(sprintf('n, edge %d (region %d)', worst_edge, crashed_region), 'Color',ink_primary,'FontWeight','normal');

nexttile;
plot(calls_sorted, R_o, 'Color', orange, 'LineWidth', 1.5); hold on;
yline(1, '--', 'Color', ink_muted);
set(gca,'Color',surface_c,'XColor',ink_muted,'YColor',ink_muted);
xlabel('RHS call index'); ylabel('R = \gamma_1\epsilon n^2 / \beta');
title('Crowding ratio R (dashed line = saturation point R=1)', 'Color',ink_primary,'FontWeight','normal');

nexttile;
plot(calls_sorted, veff_o, 'Color', aqua, 'LineWidth', 1.5);
set(gca,'Color',surface_c,'XColor',ink_muted,'YColor',ink_muted);
xlabel('RHS call index'); ylabel('v_{eff}');
title('Effective axonal velocity', 'Color',ink_primary,'FontWeight','normal');

nexttile;
plot(calls_sorted, np_o, 'Color', yellow, 'LineWidth', 1.5);
set(gca,'Color',surface_c,'XColor',ink_muted,'YColor',ink_muted);
xlabel('RHS call index'); ylabel('n''(x)');
title('RHS (nprime)', 'Color',ink_primary,'FontWeight','normal');

exportgraphics(fig, fullfile(out_dir,'axon_crash_worst_edge.png'), 'Resolution',150);
fprintf('\nSaved plot to %s\n', fullfile(out_dir,'axon_crash_worst_edge.png'));
