%% loops through ROIs, matches voxelwise neurotransmitter values to voxelwise normalized timescale values from each ROI to shared XYZ coordinates, then saves per-ROI 14-panel scatterplot PNGs (one receptor per panel) with distinct colors/markers and x-jitter for visual separation, showing only the Spearman correlation r in each panel.

% - Plain scatter (no regression line/CI)
% - Only "r=..." (Spearman) on each tile
% - Y label: "Timescale Value"
% - No title, no panel letters, no extra per-tile text
% - Per-ROI x-limits (1–99% quantiles) per receptor
% - Distinct color + marker shape per receptor
% - Density-breaking x-jitter (bigger if near-constant) so panels don't look identical
% - Save ONLY PNGs into a NEW subfolder under output_dir

%% ---- File/dir setup ----
prefix = '/mnt/scratch/NHP4CYRUS/';
output_dir = [prefix 'output/'];
if ~exist(output_dir, 'dir'); mkdir(output_dir); end
png_dir = fullfile(output_dir, 'png_scatter_corrOnly_ROIxlims_forcedDistinct'); % NEW subfolder
if ~exist(png_dir, 'dir'); mkdir(png_dir); end

%% ---- ROIs and atlases ----
rois = struct();
rois.dlPFC = struct('label', 54,  'atlas', 'C3_whole_brain.1D');
rois.latOFC= struct('label', 25,  'atlas', 'C3_whole_brain.1D');
rois.PMd   = struct('label', 81,  'atlas', 'C5_whole_brain.1D');
rois.ACC   = struct('label', 3,   'atlas', 'C3_whole_brain.1D');
rois.OFC   = struct('label', 16,  'atlas', 'C2_whole_brain.1D');
rois.LPFC  = struct('label', 50,  'atlas', 'C2_whole_brain.1D');
rois.LIP   = struct('label', 113, 'atlas', 'C4_whole_brain.1D');
rois.MT    = struct('label', 232, 'atlas', 'C4_whole_brain.1D');
rois.S2    = struct('label', 95,  'atlas', 'C4_whole_brain.1D');
rois.S1    = struct('label', 92,  'atlas', 'C3_whole_brain.1D');

%% ---- Receptor labels (TeX formatting for subscripts/Greek) ----
neuro_names_tex = { ...
    'AMPA','Kainate','NMDA', ...
    'GABA_{A}','GABA_{A}/BZ','GABA_{B}', ...
    'M_{1}','M_{2}','M_{3}', ...
    '5-HT_{1A}','5-HT_{2A}', ...
    '\alpha_{1}','\alpha_{2}','D_{1}'};

%% ---- Distinct colors (14 curated) and marker shapes (14) ----
% Colors: readable and distinct; tweak if you like.
colors = [ ...
    0.1216 0.4667 0.7059;  % blue
    1.0000 0.4980 0.0549;  % orange
    0.1725 0.6275 0.1725;  % green
    0.8392 0.1529 0.1569;  % red
    0.5804 0.4039 0.7412;  % purple
    0.5490 0.3373 0.2941;  % brown
    0.8902 0.4667 0.7608;  % pink
    0.4980 0.4980 0.4980;  % gray
    0.7373 0.7412 0.1333;  % olive
    0.0902 0.7451 0.8118;  % cyan
    0.7373 0.2235 0.3882;  % wine
    0.3176 0.7804 0.3922;  % mint
    0.9765 0.7804 0.0510;  % gold
    0.2745 0.5098 0.7059]; % steel blue
markers = {'o','^','s','d','v','>','<','p','h','o','^','s','d','v'}; % 14 entries

%% ---- Load neurotransmitter data (XYZ + 14) ----
try
    neuro_data = load([prefix 'RM_dump_nan.1D']); % 17 columns: XYZ(1-3) + 14 transmitters
catch e
    error('Failed to load RM_dump_nan.1D: %s', e.message);
end
if size(neuro_data,2) ~= 17
    error('RM_dump_nan.1D should have 17 columns, found %d.', size(neuro_data,2));
end

%% ---- Aesthetics ----
set(groot, 'defaultFigureColor','w', ...
    'defaultAxesFontName','Helvetica', ...
    'defaultAxesFontSize',10, ...
    'defaultAxesLineWidth',0.75, ...
    'defaultAxesBox','off', ...
    'defaultLineLineWidth',1.25, ...
    'defaultTextFontName','Helvetica', ...
    'defaultTextColor',[0 0 0]);

%% ---- Iterate ROIs ----
roi_names = fieldnames(rois);
for r = 1:numel(roi_names)
    roi = roi_names{r};
    label = rois.(roi).label;
    atlas_file = rois.(roi).atlas;

    % Load atlas labels
    try
        atlas = importdata(atlas_file);
    catch e
        warning('Failed to load %s for %s: %s', atlas_file, roi, e.message);
        continue;
    end
    crd_idx = find(atlas == label);
    if isempty(crd_idx)
        warning('No voxels for %s (label %d). Skipping.', roi, label);
        continue;
    end

    % Region-specific neuro data
    data_acc = neuro_data(crd_idx, :);    % XYZ + 14 NTs
    xyz_neuro_roi = data_acc(:,1:3);

    % Load and normalize timescale within ROI (winsorize 1–99%, then min-max to [0,1])
    timescale_file = [output_dir roi '_nowindows_avg.1D'];
    try
        ts_data = load(timescale_file);   % XYZ + value
    catch e
        warning('Missing timescale for %s: %s', roi, e.message);
        continue;
    end
    if size(ts_data,2) ~= 4
        warning('Timescale format unexpected for %s. Skipping.', roi);
        continue;
    end
    xyz_ts = ts_data(:,1:3);
    ts_raw = ts_data(:,4);
    lo_ts = prctile(ts_raw,1); hi_ts = prctile(ts_raw,99);
    if ~(isfinite(lo_ts) && isfinite(hi_ts)) || lo_ts == hi_ts
        lo_ts = min(ts_raw); hi_ts = max(ts_raw);
        if lo_ts == hi_ts, lo_ts = lo_ts - 0.5; hi_ts = hi_ts + 0.5; end
    end
    ts_clip = min(max(ts_raw,lo_ts),hi_ts);
    denom = max(eps, (max(ts_clip)-min(ts_clip)));
    ts_norm = (ts_clip - min(ts_clip)) / denom;   % 0..1

    % Match voxels by XYZ
    [~, idxN, idxT] = intersect(xyz_neuro_roi, xyz_ts, 'rows', 'stable');
    if isempty(idxN)
        warning('No XYZ overlap for %s. Skipping.', roi);
        continue;
    end

    X = data_acc(idxN, 4:end);   % 14 neurotransmitters (per ROI)
    Y = ts_norm(idxT);           % 0..1

    % Remove rows with any nonfinite values
    bad = any(~isfinite(X),2) | ~isfinite(Y);
    X(bad,:) = []; Y(bad) = [];
    n_vox = numel(Y);
    if n_vox < 5
        warning('Too few voxels for %s after filtering.', roi);
        continue;
    end

    % ---- Figure (no title) ----
    f = figure('Units','pixels','Position',[100 100 1300 900],'Color','w');
    t = tiledlayout(f,4,4,'TileSpacing','compact','Padding','compact'); %#ok<NASGU>

    % ---- Draw 14 small multiples with forced distinctness ----
    for i = 1:14
        ax = nexttile; %#ok<LAXES>
        xi = X(:,i);

        % ROI-specific x-limits (1–99% quantiles)
        lo = prctile(xi, 1); hi = prctile(xi, 99);
        if ~(isfinite(lo) && isfinite(hi)) || lo == hi
            lo = min(xi); hi = max(xi);
            if lo == hi, lo = lo - 0.5; hi = hi + 0.5; end
        end
        pad = 0.03*(hi - lo + eps);
        xlohi = [lo - pad, hi + pad];

        % Density-breaking x-jitter (never touches Y). 
        % Base jitter is 1% of panel width; if near-constant or ultra-dense, increase.
        panel_w = (xlohi(2) - xlohi(1));
        base_jit = 0.01 * panel_w;
        sxi = std(xi,'omitnan');
        iqrxi = iqr(xi,1);
        near_const = (sxi < max(1e-10, 1e-6*max(iqrxi, eps)));
        if near_const || n_vox > 5000
            jit = 0.035 * panel_w;  % larger separation
        elseif n_vox > 1000
            jit = 0.02  * panel_w;
        else
            jit = base_jit;
        end
        % Build display copy: tie-aware jitter—larger jitter within ties
        [~,~,ties] = unique(xi); % same value -> same tie id
        tie_counts = accumarray(ties,1);
        tie_jit = (tie_counts(ties)-1);
        xi_disp = xi + (jit .* (randn(size(xi)) .* (1 + 0.5*tie_jit)));

        % Point size/alpha by n_vox
        if n_vox > 5000
            msize = 3; malpha = 0.10;
        elseif n_vox > 1000
            msize = 4; malpha = 0.15;
        else
            msize = 8; malpha = 0.25;
        end

        % Distinct style per receptor
        col = colors(i,:);
        mk  = markers{i};

        scatter(ax, xi_disp, Y, msize, 'filled', ...
            'Marker', mk, ...
            'MarkerFaceColor', col, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', malpha);
        hold(ax, 'on');

        % Axes cosmetics
        ax.TickDir = 'out';
        ax.XLim = xlohi;
        ax.YLim = [0 1];

        % Labels
        if mod(i-1,4)==0
            yl = ylabel('Timescale Value'); 
            yl.FontSize = 9;
        end
        xl = xlabel(neuro_names_tex{i}, 'Interpreter','tex');
        xl.FontSize = 9;

        % Spearman correlation on TRUE x (no jitter)
        [rS, ~] = corr(xi, Y, 'Rows','complete', 'Type','Spearman');
        if isfinite(rS)
            text(ax, 0.98, 0.92, sprintf('r=%.2f', rS), ...
                'Units','normalized','HorizontalAlignment','right', ...
                'FontSize',9, 'BackgroundColor','w','Margin',1);
        end
    end

    drawnow;

    % ---- Save ONLY PNG into NEW subfolder ----
    base = fullfile(png_dir, ['neuro_timescale_scatterCorr_', roi]);
    try
        exportgraphics(f, [base,'.png'], 'Resolution', 600);
    catch e
        warning('exportgraphics failed for %s: %s. Using saveas fallback.', roi, e.message);
        saveas(f, [base,'.png']);
    end

    close(f);
end
