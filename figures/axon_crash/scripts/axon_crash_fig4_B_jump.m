% Figure: fsolve's shooting parameter B (edge 3) across its 8 attempts,
% showing the jump from B=+0.0088 to B=-0.092 that set n negative to begin with.
% Read-only: loads the existing debug dump, writes the figure to out_dir only.

cd('C:\Users\USER\Documents\Code\NTM_ru_ab\Transport_Network_Neumann');
out_dir = 'C:\Users\USER\Documents\Code\NTM_ru_ab\figures\axon_crash';
if ~isfolder(out_dir), mkdir(out_dir); end

d = load('axon_debug_log_region3_20261026_104500.mat');
term_log = d.term_log;
region_ids = arrayfun(@(s) s.region_i, term_log);
rows = term_log(region_ids == d.i);
calls = [rows.call];
[calls_sorted, order] = sort(calls);
rows = rows(order);

B3 = arrayfun(@(s) s.B(3), rows);
n3 = arrayfun(@(s) s.n(3), rows);

% Each attempt holds B fixed; detect attempt boundaries as changes in B.
attempt_change = [true, diff(B3) ~= 0];
attempt_start_idx = find(attempt_change);
n_attempts = numel(attempt_start_idx);
fprintf('Detected %d shooting attempts for edge 3. B per attempt:\n', n_attempts);
for a = 1:n_attempts
    fprintf('  attempt %d: B = %.6g  (starts at call %d)\n', a, B3(attempt_start_idx(a)), calls_sorted(attempt_start_idx(a)));
end

n_star = 0.00632456; % R=1 threshold, as established earlier

%% Legibility defaults
set(groot, 'DefaultAxesFontSize', 13);
set(groot, 'DefaultTextFontSize', 13);

blue   = double([0x2a 0x78 0xd6])/255;
orange = double([0xeb 0x68 0x34])/255;
ink_primary = double([0x0b 0x0b 0x0b])/255;
ink_muted   = double([0x89 0x87 0x81])/255;
surface_c   = double([0xfc 0xfc 0xfb])/255;
grid_hair   = double([0xe1 0xe0 0xd9])/255;

fig4 = figure('Color','w','Position',[100 100 1150 750]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% --- Top: B (the shooting parameter fsolve is searching over) ---
ax1 = nexttile;
stairs(calls_sorted, B3, 'Color', blue, 'LineWidth', 2.5);
hold(ax1,'on');
yline(0, '-', 'Color', ink_muted, 'LineWidth', 1, 'HandleVisibility','off');
hold(ax1,'off');
set(ax1,'Color',surface_c,'XColor',ink_muted,'YColor',ink_muted,'GridColor',grid_hair,'FontSize',13);
grid(ax1,'on');
ylabel('B (shooting parameter)','FontSize',14);
title('fsolve holds B near +0.0088, then one step overshoots past B=0', ...
    'Color',ink_primary,'FontWeight','bold','FontSize',16);
text(ax1, 0.015, 0.12, 'B = n(x=0): must be \geq 0 physically', 'Units','normalized', ...
    'Color', ink_primary, 'FontSize', 12, 'BackgroundColor',[1 1 1 0.75], 'Margin', 3);

% --- Bottom: |n| (log scale), same x-axis, threshold marked ---
ax2 = nexttile;
semilogy(calls_sorted, abs(n3), 'Color', orange, 'LineWidth', 2.5);
hold(ax2,'on');
yline(n_star, '--', 'Color', ink_muted, 'LineWidth', 1.5, 'HandleVisibility','off');
hold(ax2,'off');
set(ax2,'Color',surface_c,'XColor',ink_muted,'YColor',ink_muted,'GridColor',grid_hair,'FontSize',13);
grid(ax2,'on');
xlabel('RHS call index','FontSize',14); ylabel('|n|  (log scale)','FontSize',14);
text(ax2, 0.015, 0.92, sprintf('dashed line: n^* = %.4g (crowding saturates)', n_star), ...
    'Units','normalized', 'Color', ink_primary, 'FontSize', 12, 'BackgroundColor',[1 1 1 0.75], 'Margin', 3);

linkaxes([ax1 ax2], 'x');

exportgraphics(fig4, fullfile(out_dir,'fig4_B_jump_edge3.png'), 'Resolution',150);
fprintf('Saved fig4_B_jump_edge3.png\n');
