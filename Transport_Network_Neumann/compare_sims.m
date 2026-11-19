function compare_sims(sim_paths, labels)
% compare_sims  Compare N saved NTM_ru_ab runs: regional tau burden over time.
%   compare_sims(sim_paths)          -- labels auto-derived from filenames
%   compare_sims(sim_paths, labels)  -- labels{k} is the display name for sim_paths{k}
%
% Read-only: loads existing .mat outputs, writes figures to out_dir only.

    if nargin < 2 || isempty(labels)
        labels = cell(size(sim_paths));
        for k = 1:numel(sim_paths)
            [~, name] = fileparts(sim_paths{k});
            labels{k} = strrep(name, '_', '\_');  % escape underscore for tex interpreter
        end
    end
    assert(numel(sim_paths) == numel(labels), 'sim_paths and labels must be the same length');
    nsims = numel(sim_paths);

    out_dir = fullfile(fileparts(mfilename('fullpath')), 'figures');
    if ~isfolder(out_dir), mkdir(out_dir); end

    %% Palette (dataviz skill categorical theme, fixed order -- do not cycle)
    categorical_hex = {'#2a78d6','#eb6834','#1baf7a','#eda100','#e87ba4','#008300','#6250d6','#e34948'};
    hex2rgb = @(h) sscanf(h(2:end), '%2x%2x%2x')' / 255;
    if nsims > numel(categorical_hex)
        warning('compare_sims:tooManySeries', ...
            '%d sims requested but only %d CVD-safe categorical colors are validated; colors repeat past slot %d.', ...
            nsims, numel(categorical_hex), numel(categorical_hex));
    end
    colors = zeros(nsims,3);
    for k = 1:nsims
        colors(k,:) = hex2rgb(categorical_hex{mod(k-1,numel(categorical_hex))+1});
    end

    ink_primary   = double([0x0b 0x0b 0x0b])/255;
    ink_secondary = double([0x52 0x51 0x4e])/255;
    ink_muted     = double([0x89 0x87 0x81])/255;
    grid_hair     = double([0xe1 0xe0 0xd9])/255;
    surface_c     = double([0xfc 0xfc 0xfb])/255;

    seq_light = double([0xcd 0xe2 0xfb])/255;
    seq_dark  = double([0x0d 0x36 0x6b])/255;
    make_seq  = @(n) interp1([0 1], [seq_light; seq_dark], linspace(0,1,n));

    diverge_neg = double([0xd0 0x3b 0x3b])/255;
    diverge_pos = double([0x2a 0x78 0xd6])/255;

    %% Load all runs
    runs = cell(nsims,1);
    for k = 1:nsims
        runs{k} = load(sim_paths{k});
    end

    t = runs{1}.t;
    for k = 2:nsims
        assert(numel(t) == numel(runs{k}.t), ...
            'Run %d ("%s") has %d timesteps, run 1 has %d -- not directly comparable.', ...
            k, labels{k}, numel(runs{k}.t), numel(t));
        if ~isequal(t, runs{k}.t)
            warning('Run %d ("%s") time grid differs slightly from run 1; using run 1''s t.', k, labels{k});
        end
    end

    nregions = size(runs{1}.N, 1);
    total = cell(nsims,1);
    for k = 1:nsims
        total{k} = runs{k}.N + runs{k}.M;
        fprintf('  %-20s max(N)=%.4g max(M)=%.4g\n', labels{k}, max(runs{k}.N(:)), max(runs{k}.M(:)));
    end
    fprintf('Loaded %d runs, %d regions x %d timesteps.\n', nsims, nregions, numel(t));

    %% Figure 1: region x time heatmaps, one tile per run, shared color scale
    all_vals = cellfun(@(x) x(:), total, 'UniformOutput', false);
    clims = [min(vertcat(all_vals{:})), max(vertcat(all_vals{:}))];

    fig1 = figure('Color','w','Position',[100 100 min(420*nsims,1800) 430]);
    tiledlayout(1, nsims, 'TileSpacing','compact', 'Padding','compact');
    for k = 1:nsims
        nexttile;
        imagesc(t, 1:nregions, total{k}, clims); colormap(gca, make_seq(256));
        set(gca,'YDir','normal','Color',surface_c,'XColor',ink_muted,'YColor',ink_muted);
        xlabel('time'); if k==1, ylabel('region'); end
        title(labels{k}, 'Color', ink_primary, 'FontWeight','normal');
        cb = colorbar; cb.Label.String = 'total tau, N+M'; cb.Color = ink_secondary;
    end
    exportgraphics(fig1, fullfile(out_dir,'fig1_region_time_heatmaps.png'), 'Resolution',150);
    fprintf('Saved fig1_region_time_heatmaps.png\n');

    %% Figure 2: global tau burden over time, all N runs overlaid
    fig2 = figure('Color','w','Position',[100 100 900 550]);
    ax = axes(); hold(ax,'on');
    for k = 1:nsims
        plot(t, sum(total{k},1), 'Color', colors(k,:), 'LineWidth', 2, 'DisplayName', labels{k});
    end
    hold(ax,'off');
    set(ax,'Color',surface_c,'XColor',ink_muted,'YColor',ink_muted,'GridColor',grid_hair);
    grid(ax,'on');
    xlabel('time'); ylabel(sprintf('global total tau (summed over %d regions)', nregions));
    title('Global tau burden over time', 'Color',ink_primary,'FontWeight','normal');
    lgd = legend('Location','best'); lgd.Box = 'off'; lgd.TextColor = ink_secondary;
    exportgraphics(fig2, fullfile(out_dir,'fig2_global_burden.png'), 'Resolution',150);
    fprintf('Saved fig2_global_burden.png\n');

    %% Figure 3: per-region diff at final time, each run vs. run 1 (baseline)
    if nsims >= 2
        ncomp = nsims - 1;
        fig3 = figure('Color','w','Position',[100 100 min(500*ncomp,1800) 500]);
        tiledlayout(1, ncomp, 'TileSpacing','compact', 'Padding','compact');
        baseline_final = total{1}(:,end);
        for k = 2:nsims
            nexttile;
            final_diff = total{k}(:,end) - baseline_final;
            [sorted_diff, order] = sort(final_diff, 'descend');
            bar_colors = zeros(nregions,3);
            bar_colors(sorted_diff>=0,:) = repmat(diverge_pos, sum(sorted_diff>=0), 1);
            bar_colors(sorted_diff<0,:)  = repmat(diverge_neg, sum(sorted_diff<0), 1);
            b = bar(sorted_diff, 'FaceColor','flat');
            b.CData = bar_colors;
            set(gca,'Color',surface_c);
            title(sprintf('%s - %s', labels{k}, labels{1}), 'Color',ink_primary,'FontWeight','normal');
            xticks(1:nregions); xticklabels(string(order)); xtickangle(90);
        end
        exportgraphics(fig3, fullfile(out_dir,'fig3_per_region_diff_vs_baseline.png'), 'Resolution',150);
        fprintf('Saved fig3_per_region_diff_vs_baseline.png\n');
    end

    fprintf('\nAll figures written to %s\n', out_dir);
end