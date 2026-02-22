%% computes voxelwise timescale metrics (INT, lag-1 ACF magnitude, and first non-positive lag) in 4 segments for each monkey/area/run, averages them within each ROI for CHARM level, normalizes values per level/segment, then produces L2-aligned heatmaps (plus segment-comparison and overview plots).
function timescales_for_each_CHARM_level()
% timescales_for_each_CHARM_level
% L2-aligned intrinsic timescales plots across CHARM levels (2..6)
% - One heatmap per CHARM level per metric (int, ac_mag, lags)
% - y-axis = regions (from CHARM level 2), same order for every level
% - x-axis = 4 segments
% - higher levels fall back to L2 values if a name match can’t be found
%
% Required inputs:
%   /mnt/scratch/NHP4CYRUS/data_dump/<monkey>/<area>/ROIdump_<area><run>.1D
%   /mnt/scratch/NHP4CYRUS/CHARM_Level_2.csv ... CHARM_Level_6.csv
%   /mnt/scratch/NHP4CYRUS/1dfiles/C2_whole_brain.1D ... C6_whole_brain.1D

%% -------------------- 1) Config --------------------
prefix       = '/mnt/scratch/NHP4CYRUS';
atlas_prefix = fullfile(prefix, '1dfiles');
base_path    = fullfile(prefix, 'data_dump');

output_dir   = fullfile(prefix, 'output', 'timescales', 'comprehensive_analysis');
output_dir   = ensure_output_dir(output_dir);

no_runs   = 5;
TR        = 1.1; %#ok<NASGU>  % not used here, but you may want it later

metrics = {'int','ac_mag','lags'};   % three outputs
num_segments = 4;

%% -------------------- 2) Discover monkeys --------------------
if ~exist(base_path, 'dir')
    error('Base path does not exist: %s', base_path);
end

folder_info = dir(base_path);
folder_names = {folder_info([folder_info.isdir]).name};

% Try to exclude '.' and '..' robustly + exclude Delmar as you did
monkeys = folder_names(~ismember(folder_names, {'.','..','Delmar'}));

if isempty(monkeys)
    error('No monkey folders found under %s', base_path);
end

%% -------------------- 3) Load CHARM CSVs: region name + label per level --------------------
level_regions = repmat(struct('tbl', [], 'level', []), 1, 6);

for level = 2:6
    csv_file = fullfile(prefix, sprintf('CHARM_Level_%d.csv', level));
    if ~isfile(csv_file)
        warning('Missing %s. Skipping level %d.', csv_file, level);
        continue;
    end

    T = readtable(csv_file);

    expected_columns = {'row_id', 'name', 'level_idx'};
    if ~all(ismember(expected_columns, T.Properties.VariableNames))
        warning('CHARM_Level_%d.csv missing expected columns. Skipping.', level);
        continue;
    end

    [unique_names, ia] = unique(T.name, 'stable');

    cleaned = cell(size(unique_names));
    for k = 1:numel(unique_names)
        cleaned{k} = clean_region_name(unique_names{k});
    end

    R = table();
    R.Region = cleaned(:);
    R.Label  = T.level_idx(ia);
    R.Atlas_File = repmat({fullfile(atlas_prefix, sprintf('C%d_whole_brain.1D', level))}, height(R), 1);

    valid_mask = ~cellfun(@isempty, R.Region) & ~isnan(R.Label);
    R = R(valid_mask, :);

    level_regions(level).tbl   = R;
    level_regions(level).level = level;
end

%% -------------------- 4) Aggregate ROI values across monkeys/areas/runs --------------------
% We will build per level maps: label -> (sum over samples, count)
% for each metric, for each segment.

agg = init_agg_struct(6, metrics, num_segments);

for i = 1:numel(monkeys)
    monkey = monkeys{i};

    for area = {'cortical','subcortical'}
        area_name = area{1};

        for run = 1:no_runs
            dump_file = fullfile(base_path, monkey, area_name, ...
                sprintf('ROIdump_%s%d.1D', area_name, run));

            if ~isfile(dump_file)
                fprintf('Missing %s (skip)\n', dump_file);
                continue;
            end

            try
                data = importdata(dump_file);
            catch ME
                fprintf('Skipping %s: load failed (%s)\n', dump_file, ME.message);
                continue;
            end

            if isempty(data) || size(data,2) < 8
                % need at least 2 timepoints per segment; with 4 segments, total >= 8 is a safe minimal check
                fprintf('Skipping %s: insufficient data.\n', dump_file);
                continue;
            end

            % Compute voxel metrics for this run: (voxels x segments)
            [voxel_int, voxel_acmag, voxel_lags] = compute_timescales_4segments(data);

            % Push into each CHARM level using that level's atlas labels
            for level = 2:6
                R = level_regions(level).tbl;
                if isempty(R)
                    continue;
                end

                atlas_file = R.Atlas_File{1};
                if ~isfile(atlas_file)
                    warning('Missing atlas file %s (level %d).', atlas_file, level);
                    continue;
                end

                try
                    atlas_labels = importdata(atlas_file);
                catch
                    warning('Failed to load atlas %s (level %d).', atlas_file, level);
                    continue;
                end

                % IMPORTANT: assume atlas_labels rows correspond to ROIdump voxel rows
                num_vox = size(data,1);
                m = min(num_vox, numel(atlas_labels));
                atlas_labels = atlas_labels(1:m);

                v_int   = voxel_int(1:m, :);
                v_acmag = voxel_acmag(1:m, :);
                v_lags  = voxel_lags(1:m, :);

                % Aggregate per ROI label (only those labels present in CSV list)
                for rr = 1:height(R)
                    lbl = R.Label(rr);

                    idx = find(atlas_labels == lbl);
                    if isempty(idx)
                        continue;
                    end

                    % per-ROI mean per segment
                    roi_int   = mean(v_int(idx,:),   1, 'omitnan');
                    roi_acmag = mean(v_acmag(idx,:), 1, 'omitnan');
                    roi_lags  = mean(v_lags(idx,:),  1, 'omitnan');

                    % store into aggregator
                    agg = agg_add(agg, level, lbl, 'int',   roi_int);
                    agg = agg_add(agg, level, lbl, 'ac_mag',roi_acmag);
                    agg = agg_add(agg, level, lbl, 'lags',  roi_lags);
                end
            end
        end
    end
end

%% -------------------- 5) Build per-level tables: Region + Values --------------------
% Convert aggregated label stats into level data structs:
all_levels_data = repmat(struct('regions', [], 'values', [], 'level', []), 1, 6);

reference_regions = {};

for level = 2:6
    R = level_regions(level).tbl;
    if isempty(R)
        continue;
    end

    % extract mean values per region label from aggregator (rows = regions, cols = segments)
    [vals_int, vals_acmag, vals_lags, valid_mask] = agg_to_values(agg, level, R.Label, num_segments);

    % choose one metric at a time for heatmaps, but we keep them separately:
    level_data.int.regions   = R.Region(valid_mask);
    level_data.int.values    = vals_int(valid_mask, :);
    level_data.int.level     = level;

    level_data.ac_mag.regions = R.Region(valid_mask);
    level_data.ac_mag.values  = vals_acmag(valid_mask, :);
    level_data.ac_mag.level   = level;

    level_data.lags.regions = R.Region(valid_mask);
    level_data.lags.values  = vals_lags(valid_mask, :);
    level_data.lags.level   = level;

    % Normalize within level per segment (like your receptor normalization per column)
    level_data = normalize_level_by_column(level_data);

    % Drop all-zero rows after normalization (per metric)
    all_levels_data(level).level = level;
    all_levels_data(level).metrics = level_data;

    if level == 2
        % Use INT metric’s surviving regions as reference (you can swap to lags/ac_mag if you prefer)
        reference_regions = level_data.int.regions;
    end
end

%% -------------------- 6) Build hierarchy mapping from L2 names --------------------
if isempty(reference_regions)
    warning('No Level 2 reference regions found. Skipping aligned heatmaps.');
    fprintf('Done (no outputs).\n');
    return;
end

hierarchy_map = build_hierarchy_mapping_timescales(all_levels_data, reference_regions);

%% -------------------- 7) Create L2-aligned heatmaps per level per metric --------------------
create_hierarchical_heatmaps_timescales_L2aligned(all_levels_data, reference_regions, ...
    hierarchy_map, output_dir, metrics, num_segments);

%% -------------------- 8) Summary figures --------------------
create_segment_comparison(all_levels_data, output_dir, metrics, num_segments);
create_hierarchy_overview_timescales(all_levels_data, output_dir, metrics, num_segments);

fprintf('Timescale plots saved in "%s"\n', output_dir);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                           LOCAL FUNCTIONS                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out_dir = ensure_output_dir(out_dir)
    try
        if ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end
        tf = fullfile(out_dir, 'test_write.txt');
        fid = fopen(tf, 'w');
        if fid == -1, error('cannot write'); end
        fclose(fid);
        delete(tf);
    catch
        warning('Permission denied for %s. Using current directory instead.', out_dir);
        out_dir = fullfile(pwd, 'output', 'timescales', 'comprehensive_analysis');
        if ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end
    end
end

function cleaned_name = clean_region_name(name)
    if iscell(name), name = name{1}; end
    cleaned_name = regexprep(name, '[\\^_{}<>]|(\\[a-zA-Z]+)|(<[^>]+>)', '');
    cleaned_name = strrep(cleaned_name, '_', ' ');
    cleaned_name = strrep(cleaned_name, '-', ' ');
    cleaned_name = regexprep(cleaned_name, '([a-z])([A-Z])', '$1 $2');
    cleaned_name = regexprep(cleaned_name, '([a-zA-Z])(\d)', '$1 $2');
    cleaned_name = strtrim(regexprep(cleaned_name, '\s+', ' '));
    if ~isempty(cleaned_name)
        words = split(cleaned_name);
        for i = 1:numel(words)
            w = words{i};
            if ~isempty(w)
                words{i} = [upper(w(1)) lower(w(2:end))];
            end
        end
        cleaned_name = strjoin(words, ' ');
    end
end

function pretty = prettify_region_label(name)
    if iscell(name), name = name{1}; end
    pretty = name;
    pretty = strrep(pretty, '_', ' ');
    pretty = regexprep(pretty, '([a-z])([A-Z])', '$1 $2');
    pretty = regexprep(pretty, '([A-Z])([A-Z][a-z])', '$1 $2');
    pretty = regexprep(pretty, '([a-zA-Z])(\d)', '$1 $2');
    pretty = strtrim(regexprep(pretty, '\s+', ' '));
end

function [v_int, v_acmag, v_lags] = compute_timescales_4segments(data)
    % data: voxels x time
    num_voxels = size(data,1);
    T = size(data,2);

    window_size = floor(T/4);
    v_int   = NaN(num_voxels, 4);
    v_acmag = NaN(num_voxels, 4);
    v_lags  = NaN(num_voxels, 4);

    for win = 1:4
        start_idx = (win-1)*window_size + 1;
        end_idx   = start_idx + window_size - 1;
        if win == 4
            end_idx = T; % include remainder in last window
        end
        segment = data(:, start_idx:end_idx);

        if size(segment,2) < 2
            continue;
        end

        max_lag = size(segment,2) - 1;

        for voxel = 1:num_voxels
            try
                acf = autocorr(segment(voxel,:), max_lag);

                % int: sum positive acf (including lag0)
                pos = acf > 0;
                v_int(voxel, win) = sum(acf(pos));

                % ac_mag: lag1
                if numel(acf) >= 2
                    v_acmag(voxel, win) = acf(2);
                else
                    v_acmag(voxel, win) = NaN;
                end

                % lags: first non-positive (<=0), last positive lag index
                first_nonpos = find(~pos, 1, 'first');
                if isempty(first_nonpos)
                    v_lags(voxel, win) = max_lag;
                else
                    v_lags(voxel, win) = first_nonpos - 1;
                end

            catch
                % leave NaNs
                continue;
            end
        end
    end
end

function agg = init_agg_struct(max_level, metrics, num_segments)
    agg = repmat(struct(), 1, max_level);
    for level = 1:max_level
        for m = 1:numel(metrics)
            metric = metrics{m};
            agg(level).(metric) = containers.Map('KeyType','double','ValueType','any');
        end
        agg(level).num_segments = num_segments;
    end
end

function agg = agg_add(agg, level, label, metric, roi_vec)
    mp = agg(level).(metric);

    if ~isKey(mp, label)
        mp(label) = struct('sum', zeros(1, numel(roi_vec)), 'count', zeros(1, numel(roi_vec)));
    end

    st = mp(label);

    good = ~isnan(roi_vec);
    st.sum(good)   = st.sum(good)   + roi_vec(good);
    st.count(good) = st.count(good) + 1;

    mp(label) = st;
    agg(level).(metric) = mp;
end

function [vals_int, vals_acmag, vals_lags, valid_mask] = agg_to_values(agg, level, labels, num_segments)
    vals_int   = NaN(numel(labels), num_segments);
    vals_acmag = NaN(numel(labels), num_segments);
    vals_lags  = NaN(numel(labels), num_segments);

    valid_mask = false(numel(labels), 1);

    for i = 1:numel(labels)
        lbl = labels(i);

        % int
        if isKey(agg(level).int, lbl)
            st = agg(level).int(lbl);
            v = st.sum ./ max(st.count, 1);
            vals_int(i,:) = v;
        end

        % ac_mag
        if isKey(agg(level).ac_mag, lbl)
            st = agg(level).ac_mag(lbl);
            v = st.sum ./ max(st.count, 1);
            vals_acmag(i,:) = v;
        end

        % lags
        if isKey(agg(level).lags, lbl)
            st = agg(level).lags(lbl);
            v = st.sum ./ max(st.count, 1);
            vals_lags(i,:) = v;
        end

        if any(~isnan(vals_int(i,:))) || any(~isnan(vals_acmag(i,:))) || any(~isnan(vals_lags(i,:)))
            valid_mask(i) = true;
        end
    end
end

function level_data = normalize_level_by_column(level_data)
    % Normalize each metric separately, per segment column, to [0,1]
    metric_list = fieldnames(level_data);
    for mm = 1:numel(metric_list)
        metric = metric_list{mm};
        if ~isstruct(level_data.(metric)), continue; end

        V = level_data.(metric).values;
        for c = 1:size(V,2)
            col = V(:,c);
            good = ~isnan(col);
            if any(good)
                mn = min(col(good));
                mx = max(col(good));
                if mn ~= mx
                    col(good) = (col(good) - mn) ./ (mx - mn);
                else
                    col(good) = 0;
                end
                col(~good) = 0;
            else
                col(:) = 0;
            end
            V(:,c) = col;
        end

        % drop all-zero rows
        keep = sum(V,2) > 0;
        level_data.(metric).regions = level_data.(metric).regions(keep);
        level_data.(metric).values  = V(keep,:);
    end
end

function hierarchy_map = build_hierarchy_mapping_timescales(all_levels_data, reference_regions)
    nL2 = numel(reference_regions);
    hierarchy_map = struct([]);
    for i = 1:nL2
        ref_region = reference_regions{i};
        hierarchy_map(i).level2_region = ref_region;
        hierarchy_map(i).children = cell(5, 1); % slots for levels 2..6
        hierarchy_map(i).children{1} = {ref_region};

        ref_clean = lower(strrep(ref_region, ' ', ''));

        for level = 3:6
            slot = level - 1; % 3->2 ... 6->5
            if numel(all_levels_data) < level || isempty(all_levels_data(level).metrics)
                hierarchy_map(i).children{slot} = {};
                continue;
            end

            % use int metric’s regions list as the “available regions” set for matching
            if ~isfield(all_levels_data(level).metrics, 'int') || isempty(all_levels_data(level).metrics.int.regions)
                hierarchy_map(i).children{slot} = {};
                continue;
            end

            level_regions = all_levels_data(level).metrics.int.regions;
            matches = {};
            for j = 1:numel(level_regions)
                child_region = level_regions{j};
                child_clean  = lower(strrep(child_region, ' ', ''));
                if contains(child_clean, ref_clean) || contains(ref_clean, child_clean)
                    matches{end+1} = child_region; %#ok<AGROW>
                end
            end

            if isempty(matches)
                if any(strcmp(level_regions, ref_region))
                    matches = {ref_region};
                end
            end

            hierarchy_map(i).children{slot} = matches;
        end
    end
end

function idx = find_best_region_match(target_name, level_regions)
    idx = [];
    if isempty(level_regions), return; end

    tgt = lower(regexprep(target_name, '[^a-z0-9]', ''));
    best_score = 0;

    for k = 1:numel(level_regions)
        cand = level_regions{k};
        cand_norm = lower(regexprep(cand, '[^a-z0-9]', ''));

        if strcmp(tgt, cand_norm)
            idx = k; return;
        end

        if contains(cand_norm, tgt) || contains(tgt, cand_norm)
            score = min(numel(cand_norm), numel(tgt));
            if score > best_score
                best_score = score;
                idx = k;
            end
        end
    end
end

function create_hierarchical_heatmaps_timescales_L2aligned(all_levels_data, reference_regions, hierarchy_map, output_dir, metrics, num_segments)
    seg_labels = arrayfun(@(k) sprintf('Segment %d', k), 1:num_segments, 'UniformOutput', false);
    num_L2 = numel(reference_regions);

    for m = 1:numel(metrics)
        metric = metrics{m};

        % L2 fallback struct
        L2_vals = [];
        if numel(all_levels_data) >= 2 && ~isempty(all_levels_data(2).metrics) && isfield(all_levels_data(2).metrics, metric)
            L2_vals = all_levels_data(2).metrics.(metric);
        end

        for level = 2:6
            if numel(all_levels_data) < level || isempty(all_levels_data(level).metrics) || ~isfield(all_levels_data(level).metrics, metric)
                continue;
            end

            if isempty(all_levels_data(level).metrics.(metric).regions)
                continue;
            end

            slot = level - 1;

            aligned = zeros(num_L2, num_segments);
            display_regions = reference_regions;

            level_regions = all_levels_data(level).metrics.(metric).regions;
            level_values  = all_levels_data(level).metrics.(metric).values;

            for i = 1:num_L2
                l2_name = reference_regions{i};
                children = hierarchy_map(i).children{slot};

                row_val = [];

                % 1) explicit children
                if ~isempty(children)
                    child_vals = [];
                    for j = 1:numel(children)
                        child_region = children{j};
                        idx = find(strcmp(level_regions, child_region));
                        if ~isempty(idx)
                            child_vals = [child_vals; level_values(idx,:)]; %#ok<AGROW>
                        end
                    end
                    if ~isempty(child_vals)
                        row_val = mean(child_vals, 1);
                    end
                end

                % 2) fuzzy match
                if isempty(row_val)
                    best_idx = find_best_region_match(l2_name, level_regions);
                    if ~isempty(best_idx)
                        row_val = level_values(best_idx,:);
                    end
                end

                % 3) fallback to L2
                if isempty(row_val) && ~isempty(L2_vals)
                    l2_idx = find(strcmp(L2_vals.regions, l2_name), 1);
                    if ~isempty(l2_idx)
                        row_val = L2_vals.values(l2_idx,:);
                    end
                end

                % 4) zeros
                if isempty(row_val)
                    row_val = zeros(1, num_segments);
                end

                aligned(i,:) = row_val;
            end

            % prettify y labels
            ylab = cell(size(display_regions));
            for i = 1:numel(display_regions)
                pr = prettify_region_label(display_regions{i});
                if numel(pr) > 50, pr = [pr(1:47) '...']; end
                ylab{i} = pr;
            end

            fig_height = max(600, min(1400, 35*num_L2));
            f = figure('Position', [50 50 1200 fig_height], 'Color','w');

            imagesc(aligned);
            colormap(parula);
            cb = colorbar;
            cb.Label.String = sprintf('Normalized %s', metric);
            cb.FontSize = 12;
            cb.Label.FontWeight = 'bold';

            xticks(1:num_segments);
            xticklabels(seg_labels);
            xlabel('Segment', 'FontSize', 13, 'FontWeight','bold');

            yticks(1:numel(ylab));
            yticklabels(ylab);
            ylabel('Regions (Level 2 reference)', 'FontSize', 13, 'FontWeight','bold');

            title(sprintf('CHARM Level %d: Timescales (%s) (L2-aligned)', level, metric), ...
                'FontSize', 14, 'FontWeight','bold');

            set(gca, 'FontSize', 10, 'LineWidth', 1, 'TickLabelInterpreter','tex');
            axis tight;

            filename = sprintf('Level_%d_%s_Heatmap', level, metric);
            png_file = fullfile(output_dir, [filename '.png']);
            pdf_file = fullfile(output_dir, [filename '.pdf']);

            try
                exportgraphics(f, png_file, 'Resolution', 300);
                exportgraphics(f, pdf_file, 'ContentType', 'vector');
            catch
                saveas(f, png_file);
            end

            close(f);
        end
    end
end

function create_segment_comparison(all_levels_data, output_dir, metrics, num_segments)
    seg_labels = arrayfun(@(k) sprintf('S%d', k), 1:num_segments, 'UniformOutput', false);

    for m = 1:numel(metrics)
        metric = metrics{m};

        % collect per level segment means
        levels = 2:6;
        M = NaN(numel(levels), num_segments);

        for i = 1:numel(levels)
            level = levels(i);
            if numel(all_levels_data) < level || isempty(all_levels_data(level).metrics) || ~isfield(all_levels_data(level).metrics, metric)
                continue;
            end
            V = all_levels_data(level).metrics.(metric).values;
            if isempty(V), continue; end
            for s = 1:num_segments
                col = V(:,s);
                col = col(~isnan(col));
                if isempty(col), M(i,s) = 0; else, M(i,s) = mean(col); end
            end
        end

        f = figure('Position', [50 50 1200 700], 'Color','w');
        bar(M, 'grouped');
        xticks(1:numel(levels));
        xticklabels(arrayfun(@(x) sprintf('L%d', x), levels, 'UniformOutput', false));
        xlabel('CHARM Level', 'FontSize', 12, 'FontWeight','bold');
        ylabel(sprintf('Mean normalized %s', metric), 'FontSize', 12, 'FontWeight','bold');
        title(sprintf('Segment Comparison Across CHARM Levels (%s)', metric), 'FontSize', 14, 'FontWeight','bold');
        legend(seg_labels, 'Location','best');
        grid on;

        png_file = fullfile(output_dir, sprintf('Segment_Comparison_%s.png', metric));
        pdf_file = fullfile(output_dir, sprintf('Segment_Comparison_%s.pdf', metric));
        try
            exportgraphics(f, png_file, 'Resolution', 300);
            exportgraphics(f, pdf_file, 'ContentType','vector');
        catch
            saveas(f, png_file);
        end
        close(f);
    end
end

function create_hierarchy_overview_timescales(all_levels_data, output_dir, metrics, num_segments)
    levels = 2:6;

    for m = 1:numel(metrics)
        metric = metrics{m};

        % matrix: rows = levels, cols = segments
        H = NaN(numel(levels), num_segments);

        for i = 1:numel(levels)
            level = levels(i);
            if numel(all_levels_data) < level || isempty(all_levels_data(level).metrics) || ~isfield(all_levels_data(level).metrics, metric)
                continue;
            end
            V = all_levels_data(level).metrics.(metric).values;
            if isempty(V), continue; end
            for s = 1:num_segments
                col = V(:,s);
                col = col(~isnan(col));
                if isempty(col), H(i,s) = 0; else, H(i,s) = mean(col); end
            end
        end
        H(isnan(H)) = 0;

        f = figure('Position', [50 50 900 600], 'Color','w');
        imagesc(H);
        colormap(parula);
        cb = colorbar;
        cb.Label.String = sprintf('Mean normalized %s', metric);

        xticks(1:num_segments);
        xticklabels(arrayfun(@(k) sprintf('Segment %d', k), 1:num_segments, 'UniformOutput', false));
        yticks(1:numel(levels));
        yticklabels(arrayfun(@(x) sprintf('Level %d', x), levels, 'UniformOutput', false));

        xlabel('Segment', 'FontSize', 12, 'FontWeight','bold');
        ylabel('CHARM Level', 'FontSize', 12, 'FontWeight','bold');
        title(sprintf('Hierarchy Overview: %s Across Levels (Segment Means)', metric), 'FontSize', 14, 'FontWeight','bold');

        axis tight;

        png_file = fullfile(output_dir, sprintf('Hierarchy_Overview_%s.png', metric));
        pdf_file = fullfile(output_dir, sprintf('Hierarchy_Overview_%s.pdf', metric));
        try
            exportgraphics(f, png_file, 'Resolution', 300);
            exportgraphics(f, pdf_file, 'ContentType','vector');
        catch
            saveas(f, png_file);
        end
        close(f);
    end
end
