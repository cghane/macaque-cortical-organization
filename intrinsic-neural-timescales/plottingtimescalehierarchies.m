%% plotting_timescales_hierarchical_norm01.m
% Hierarchical representation (one plot per segment, regions on x)
% Timescales are MIN–MAX normalized to [0,1] within each SEGMENT across regions.

%% =========================
%  Groups
%  =========================
group1 = {'S1', 'S2'};
group2 = {'MT', 'LIP', 'OFC', 'LPFC', 'ACC'};
group3 = {'PMd', 'latOFC', 'dlPFC'};

% Color schemes for each group (used for per-region markers)
colors1 = {[0.5, 0, 0.5], [0.7, 0, 0.7]};  % Purples for S1-S2
colors2 = {[0, 0.7, 0], [0, 0.6, 0.2], [0, 0.5, 0.4], [0, 0.4, 0.6], [0, 0.3, 0.8]};  % Greens for MT-ACC
colors3 = {[0.8, 0.4, 0], [0.8, 0.2, 0.2], [0.8, 0, 0.4]};  % Oranges/Reds for PMd group

% Robust output directory (handles name collisions & permissions)
outdir = ensure_outdir('output');

% Global visual style
set(0,'DefaultFigureColor','w', ...
      'DefaultAxesFontName','Arial', ...
      'DefaultAxesFontSize',10, ...
      'DefaultAxesLineWidth',0.75, ...
      'DefaultLineLineWidth',1.25);

%% =========================
%  Generate plots (normalized 0–1, pretty)
%  =========================
plot_group_norm01(group1, colors1, 'S1_S2', outdir);
plot_group_norm01(group2, colors2, 'MT_ACC', outdir);
plot_group_norm01(group3, colors3, 'PMd_dlPFC', outdir);

fprintf('Plots (normalized 0–1) saved in "%s" as PNG/PDF/FIG.\n', outdir);

%% =======================================================================
%                               FUNCTIONS
% ========================================================================

function plot_group_norm01(group, colors, group_name, outdir)
    nSeg = 4;
    nReg = numel(group);
    raw = nan(nReg, nSeg);

    % -------- load raw values for this group (col 4 mean, as before) -------
    for seg = 1:nSeg
        for r = 1:nReg
            region = group{r};
            try
                S = load(sprintf('%s_avg_segment%d.mat', region, seg));
                raw(r, seg) = mean(S.avgMatrix(:, 4), 'omitnan');
            catch ME
                fprintf('Error loading %s segment %d: %s\n', region, seg, ME.message);
            end
        end
    end

    % -------- FIXED: min–max normalize within each SEGMENT across regions ---------
    norm01 = nan(size(raw));
    mn = min(raw, [], 1, 'omitnan');   % per-segment min across regions
    mx = max(raw, [], 1, 'omitnan');   % per-segment max across regions
    span = mx - mn;

    for seg = 1:nSeg
        if ~isfinite(span(seg)) || span(seg) == 0
            norm01(:, seg) = 0.5;
        else
            norm01(:, seg) = (raw(:, seg) - mn(seg)) ./ span(seg);
        end
    end

    % -------- plotting per segment (hierarchical layout preserved) ---------
    ylims = [0 1];

    for seg = 1:nSeg
        y = norm01(:, seg);

        f = figure('Position',[100,100,820,580],'Renderer','painters'); hold on;

        plot(1:nReg, y, '-', 'Color',[0.6 0.6 0.6], 'LineWidth',1);

        for r = 1:nReg
            mc = colors{min(r, numel(colors))};
            plot(r, y(r), 'o', ...
                'MarkerFaceColor', mc, ...
                'MarkerEdgeColor', 'k', ...
                'MarkerSize', 7, ...
                'LineWidth', 0.75);
        end

        xlim([0.5, nReg+0.5]); ylim(ylims);
        grid on; ax = gca;
        ax.GridAlpha = 0.15; ax.MinorGridAlpha = 0.07; ax.YGrid = 'on'; ax.XGrid = 'off';
        xticks(1:nReg); xticklabels(group); xtickangle(25);
        ylabel('Normalized timescale (0–1)'); xlabel('Region');
        title(sprintf('%s — Segment %d', group_name, seg), 'FontWeight','bold');

        set(gca,'Box','off','TickDir','out','TickLength',[0.010 0.010]);

        base = fullfile(outdir, sprintf('%s_segment%d_norm01', group_name, seg));
        export_block(f, base);
        close(f);
    end
end

function outdir = ensure_outdir(baseName)
    if exist(baseName, 'file') && ~exist(baseName, 'dir')
        warning('"%s" exists as a file. Creating "%s_figs" instead.', baseName, [baseName '_figs']);
        baseName = [baseName '_figs'];
    end
    if ~exist(baseName, 'dir')
        mkdir(baseName);
    end
    testFile = fullfile(baseName, sprintf('._writetest_%s.txt', char(java.util.UUID.randomUUID)));
    fid = fopen(testFile, 'w');
    if fid == -1
        error('Cannot write to "%s". Check permissions.', baseName);
    else
        fclose(fid); delete(testFile);
    end
    outdir = baseName;
end

function export_block(figHandle, basepath)
    pngPath = [basepath '.png'];
    pdfPath = [basepath '.pdf'];
    figPath = [basepath '.fig'];

    try exportgraphics(figHandle, pngPath, 'Resolution', 300); catch, print(figHandle, pngPath, '-dpng', '-r300'); end
    try exportgraphics(figHandle, pdfPath,  'ContentType', 'vector'); catch, print(figHandle, pdfPath, '-dpdf', '-painters'); end
    try savefig(figHandle, figPath); catch, end
end
