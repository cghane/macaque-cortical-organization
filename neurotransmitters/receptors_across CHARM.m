%% computes ROI-mean receptor densities for every region at CHARM levels 2–6, min–max normalizes each receptor column within each level, then builds L2-aligned receptor heatmaps (with fuzzy name matching and L2 fallback) plus system-level and hierarchy-overview summary plots.
function receptors_for_each_CHARM_level()
% receptors_for_each_CHARM_level
% L2-aligned neurotransmitter receptor plots across CHARM levels (2..6)
% - one heatmap per CHARM level
% - y-axis = regions (from CHARM level 2), same order for every level
% - higher levels fall back to L2 values if a name match can’t be found
% - region labels prettified
%
% Make sure these exist:
%   /mnt/scratch/NHP4CYRUS/RM_dump_nan.1D
%   /mnt/scratch/NHP4CYRUS/CHARM_Level_2.csv ... CHARM_Level_6.csv
%   /mnt/scratch/NHP4CYRUS/1dfiles/C2_whole_brain.1D ... C6_whole_brain.1D

%% 1. Define receptor names
receptor_names = { ...
    'AMPA_Glutamate', 'Kainate_Glutamate', 'NMDA_Glutamate', ...
    'GABA_A', 'GABA_A_BZ', 'GABA_B', ...
    'M1_Acetylcholine', 'M2_Acetylcholine', 'M3_Acetylcholine', ...
    'HT1A_Serotonin', 'HT2A_Serotonin', ...
    'Alpha1_Noradrenaline', 'Alpha2_Noradrenaline', ...
    'D1_Dopamine'};

%% 2. Paths
prefix       = '/mnt/scratch/NHP4CYRUS';
atlas_prefix = '/mnt/scratch/NHP4CYRUS/1dfiles';
output_dir   = fullfile(prefix, 'output', 'neurotransmitters', 'comprehensive_analysis');

%% 3. Check / create output dir
output_dir = ensure_output_dir(output_dir);

%% 4. Load receptor matrix
rm_file = fullfile(prefix, 'RM_dump_nan.1D');
if ~isfile(rm_file)
    error('File %s does not exist or is not accessible.', rm_file);
end
neuro_data = load(rm_file);
if size(neuro_data, 2) ~= 17
    error('RM_dump_nan.1D should have 17 columns, found %d.', size(neuro_data, 2));
end

%% 5. Prepare storage
all_levels_data = repmat(struct('regions', [], 'values', [], 'level', []), 1, 6);
reference_regions = {};

%% 6. Process CHARM levels 2..6
for level = 2:6
    csv_file = fullfile(prefix, sprintf('CHARM_Level_%d.csv', level));
    if ~isfile(csv_file)
        warning('Missing %s. Skipping level %d.', csv_file, level);
        continue;
    end

    try
        level_table = readtable(csv_file);
    catch ME
        warning('Failed to read %s: %s. Skipping level %d.', csv_file, ME.message, level);
        continue;
    end

    expected_columns = {'row_id', 'name', 'level_idx'};
    if ~all(ismember(expected_columns, level_table.Properties.VariableNames))
        warning('CHARM_Level_%d.csv does not have expected columns. Skipping.', level);
        continue;
    end

    % unique names for that level
    [unique_names, ia] = unique(level_table.name, 'stable');

    level_regions = table();
    cleaned_names = cell(size(unique_names));
    for k = 1:numel(unique_names)
        cleaned_names{k} = clean_region_name(unique_names{k});
    end
    level_regions.Region = cleaned_names(:);
    level_regions.Label  = level_table.level_idx(ia);

    % attach atlas file for this level
    atlas_file = fullfile(atlas_prefix, sprintf('C%d_whole_brain.1D', level));
    level_regions.Atlas_File = repmat({atlas_file}, height(level_regions), 1);

    % keep valid
    valid_mask = ~cellfun(@isempty, level_regions.Region) & ~isnan(level_regions.Label);
    level_regions = level_regions(valid_mask, :);

    if isempty(level_regions)
        warning('No valid regions for CHARM level %d. Skipping.', level);
        continue;
    end

    % compute means
    [mean_values, valid_regions] = compute_mean_values(level_regions, neuro_data, receptor_names);

    if isempty(mean_values)
        warning('No usable data for CHARM level %d. Skipping.', level);
        continue;
    end

    % normalize per receptor column
    for r = 1:numel(receptor_names)
        col = mean_values(:, r);
        good = ~isnan(col);
        if any(good)
            mn = min(col(good));
            mx = max(col(good));
            if mn ~= mx
                col(good) = (col(good) - mn) ./ (mx - mn);
            else
                col(good) = 0;
            end
            mean_values(:, r) = col;
        else
            mean_values(:, r) = 0;
        end
    end

    % drop all-zero rows
    row_sum = sum(mean_values, 2);
    keep_rows = row_sum > 0;
    if ~any(keep_rows)
        warning('All regions are zero after normalization for level %d. Skipping.', level);
        continue;
    end

    filtered_regions = valid_regions.Region(keep_rows);
    filtered_values  = mean_values(keep_rows, :);
    filtered_values(isnan(filtered_values)) = 0;

    % store
    all_levels_data(level).regions = filtered_regions;
    all_levels_data(level).values  = filtered_values;
    all_levels_data(level).level   = level;

    % level 2 = reference
    if level == 2
        reference_regions = filtered_regions;
    end
end

%% 7. Build hierarchical / aligned heatmaps
if ~isempty(reference_regions)
    hierarchy_map = build_hierarchy_mapping(all_levels_data, reference_regions);
    create_hierarchical_heatmaps_L2aligned(all_levels_data, reference_regions, ...
        hierarchy_map, receptor_names, output_dir);
else
    warning('No Level 2 reference regions were found. Skipping aligned heatmaps.');
end

%% 8. Summary figures
create_system_comparison(all_levels_data, receptor_names, output_dir);
create_hierarchy_overview(all_levels_data, receptor_names, output_dir);

fprintf('Plots saved in "%s"\n', output_dir);

end  % <-- end main function


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
        if fid == -1
            error('cannot write');
        end
        fclose(fid);
        delete(tf);
    catch
        warning('Permission denied for %s. Using current directory instead.', out_dir);
        out_dir = fullfile(pwd, 'output', 'neurotransmitters', 'comprehensive_analysis');
        if ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end
    end
end


function hierarchy_map = build_hierarchy_mapping(all_levels_data, reference_regions)
    % For each L2 region, try to find similarly named regions in levels 3..6
    nL2 = numel(reference_regions);
    hierarchy_map = struct([]);
    for i = 1:nL2
        ref_region = reference_regions{i};
        hierarchy_map(i).level2_region = ref_region;
        hierarchy_map(i).children = cell(5, 1);   % slots for levels 2..6
        hierarchy_map(i).children{1} = {ref_region};

        ref_clean = lower(strrep(ref_region, ' ', ''));

        for level = 3:6
            slot = level - 1;   % level 3 -> 2, ..., level 6 -> 5
            if numel(all_levels_data) < level
                hierarchy_map(i).children{slot} = {};
                continue;
            end
            if isempty(all_levels_data(level).regions)
                hierarchy_map(i).children{slot} = {};
                continue;
            end

            level_regions = all_levels_data(level).regions;
            matches = {};
            for j = 1:numel(level_regions)
                child_region = level_regions{j};
                child_clean  = lower(strrep(child_region, ' ', ''));
                if contains(child_clean, ref_clean) || contains(ref_clean, child_clean)
                    matches{end+1} = child_region; %#ok<AGROW>
                end
            end

            if isempty(matches)
                % exact name present?
                if any(strcmp(level_regions, ref_region))
                    matches = {ref_region};
                end
            end

            hierarchy_map(i).children{slot} = matches;
        end
    end
end


function create_hierarchical_heatmaps_L2aligned(all_levels_data, reference_regions, hierarchy_map, receptor_names, output_dir)
    % One figure per CHARM level, all with L2 regions on y-axis
    display_receptor_names = { ...
        'AMPA (Glutamate)', 'Kainate (Glutamate)', 'NMDA (Glutamate)', ...
        'GABA_A', 'GABA_{ABZ}', 'GABA_B', ...
        'M_1 (ACh)', 'M_2 (ACh)', 'M_3 (ACh)', ...
        '5-HT_{1A} (Serotonin)', '5-HT_{2A} (Serotonin)', ...
        '\alpha_1 (Noradrenaline)', '\alpha_2 (Noradrenaline)', ...
        'D_1 (Dopamine)'};

    num_receptors  = numel(receptor_names);
    num_L2_regions = numel(reference_regions);

    % L2 fallback
    L2_vals = [];
    if numel(all_levels_data) >= 2
        if ~isempty(all_levels_data(2).regions)
            L2_vals = all_levels_data(2);
        end
    end

    for level = 2:6
        if numel(all_levels_data) < level
            continue;
        end
        if isempty(all_levels_data(level).regions)
            continue;
        end

        slot = level - 1;
        aligned_values  = zeros(num_L2_regions, num_receptors);
        display_regions = reference_regions;

        for i = 1:num_L2_regions
            l2_name = reference_regions{i};
            children = hierarchy_map(i).children{slot};

            row_val = [];

            % 1) try explicit children
            if ~isempty(children)
                child_values = [];
                for j = 1:numel(children)
                    child_region = children{j};
                    idx = find(strcmp(all_levels_data(level).regions, child_region));
                    if ~isempty(idx)
                        child_values = [child_values; all_levels_data(level).values(idx, :)]; %#ok<AGROW>
                    end
                end
                if ~isempty(child_values)
                    row_val = mean(child_values, 1);
                end
            end

            % 2) fuzzy match inside this level
            if isempty(row_val)
                best_idx = find_best_region_match(l2_name, all_levels_data(level).regions);
                if ~isempty(best_idx)
                    row_val = all_levels_data(level).values(best_idx, :);
                end
            end

            % 3) fallback to L2 value
            if isempty(row_val)
                if ~isempty(L2_vals)
                    l2_idx = find(strcmp(L2_vals.regions, l2_name), 1);
                    if ~isempty(l2_idx)
                        row_val = L2_vals.values(l2_idx, :);
                    end
                end
            end

            % 4) if still nothing, use zeros
            if isempty(row_val)
                row_val = zeros(1, num_receptors);
            end

            aligned_values(i, :) = row_val;
        end

        % prettify y labels
        display_regions_short = cell(size(display_regions));
        for i = 1:numel(display_regions)
            pr = prettify_region_label(display_regions{i});
            if numel(pr) > 50
                pr = [pr(1:47) '...'];
            end
            display_regions_short{i} = pr;
        end

        % figure height
        fig_height = max(600, min(1400, 35 * num_L2_regions));
        f = figure('Position', [50 50 1400 fig_height], 'Color', 'w');

        imagesc(aligned_values);
        colormap(parula);
        cb = colorbar;
        cb.Label.String = 'Normalized Receptor Density';
        cb.FontSize = 12;
        cb.Label.FontWeight = 'bold';

        xticks(1:numel(display_receptor_names));
        xticklabels(display_receptor_names);
        xtickangle(45);
        xlabel('Receptor Type', 'FontSize', 13, 'FontWeight', 'bold');

        yticks(1:numel(display_regions_short));
        yticklabels(display_regions_short);
        ylabel('Regions', 'FontSize', 13, 'FontWeight', 'bold');

        title(sprintf('CHARM Level %d: Neurotransmitter Receptor Distribution (L2-aligned)', level), ...
            'FontSize', 14, 'FontWeight', 'bold');

        set(gca, 'FontSize', 10, 'LineWidth', 1, 'TickLabelInterpreter', 'tex');
        axis tight;

        filename = sprintf('Level_%d_Heatmap', level);
        png_file  = fullfile(output_dir, [filename '.png']);
        pdf_file  = fullfile(output_dir, [filename '.pdf']);
        try
            exportgraphics(f, png_file, 'Resolution', 300);
            exportgraphics(f, pdf_file, 'ContentType', 'vector');
        catch
            saveas(f, png_file);
        end

        close(f);
    end
end


function idx = find_best_region_match(target_name, level_regions)
    % very simple fuzzy matcher
    idx = [];
    if isempty(level_regions)
        return;
    end

    tgt = lower(regexprep(target_name, '[^a-z0-9]', ''));
    best_score = 0;

    for k = 1:numel(level_regions)
        cand = level_regions{k};
        cand_norm = lower(regexprep(cand, '[^a-z0-9]', ''));

        if strcmp(tgt, cand_norm)
            idx = k;
            return;
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


function cleaned_name = clean_region_name(name)
    if iscell(name)
        name = name{1};
    end

    % remove LaTeX-ish chars
    cleaned_name = regexprep(name, '[\\^_{}<>]|(\\[a-zA-Z]+)|(<[^>]+>)', '');

    % underscores, hyphens -> space
    cleaned_name = strrep(cleaned_name, '_', ' ');
    cleaned_name = strrep(cleaned_name, '-', ' ');

    % camelCase split
    cleaned_name = regexprep(cleaned_name, '([a-z])([A-Z])', '$1 $2');

    % letter-number split
    cleaned_name = regexprep(cleaned_name, '([a-zA-Z])(\d)', '$1 $2');

    % collapse spaces
    cleaned_name = strtrim(regexprep(cleaned_name, '\s+', ' '));

    % title case
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
    if iscell(name)
        name = name{1};
    end

    pretty = name;

    % underscores -> space
    pretty = strrep(pretty, '_', ' ');

    % split camel / Pascal
    pretty = regexprep(pretty, '([a-z])([A-Z])', '$1 $2');
    pretty = regexprep(pretty, '([A-Z])([A-Z][a-z])', '$1 $2');

    % letters followed by digit
    pretty = regexprep(pretty, '([a-zA-Z])(\d)', '$1 $2');

    % collapse multiple spaces
    pretty = strtrim(regexprep(pretty, '\s+', ' '));
end


function [mean_values, valid_regions] = compute_mean_values(level_regions, neuro_data, receptor_names)
    num_receptors = numel(receptor_names);
    n_regions     = height(level_regions);

    mean_values   = NaN(n_regions, num_receptors);
    valid_regions = table();
    valid_regions.Region = {};
    valid_mask    = false(n_regions, 1);

    for r = 1:n_regions
        label     = level_regions.Label(r);
        atlasfile = level_regions.Atlas_File{r};

        if ~isfile(atlasfile)
            continue;
        end

        try
            atlas_data = importdata(atlasfile);
        catch
            continue;
        end

        voxel_idx = find(atlas_data == label);
        if isempty(voxel_idx)
            continue;
        end

        data_acc = neuro_data(voxel_idx, :);

        for k = 1:num_receptors
            col_idx = k + 3;  % receptor columns start at col 4
            if size(data_acc, 2) >= col_idx
                mean_values(r, k) = mean(data_acc(:, col_idx), 'omitnan');
            end
        end

        valid_mask(r) = true;
        valid_regions.Region{end+1, 1} = level_regions.Region{r}; %#ok<AGROW>
    end

    mean_values   = mean_values(valid_mask, :);
    valid_regions = valid_regions(~cellfun(@isempty, valid_regions.Region), :);
end


function create_system_comparison(all_levels_data, receptor_names, output_dir)
    systems = {'Glutamate', 'GABA', 'Acetylcholine', 'Serotonin', 'Noradrenaline', 'Dopamine'};

    f = figure('Position', [50 50 1600 1000], 'Color', 'w');

    for s = 1:numel(systems)
        subplot(2, 3, s);
        sys = systems{s};

        sys_idx = contains(receptor_names, sys);

        level_means  = [];
        level_labels = {};

        for level = 2:6
            if numel(all_levels_data) < level
                continue;
            end
            if isempty(all_levels_data(level).values)
                continue;
            end
            vals = all_levels_data(level).values(:, sys_idx);
            vals = vals(~isnan(vals));
            if isempty(vals)
                level_means(end+1, 1) = 0; %#ok<AGROW>
            else
                level_means(end+1, 1) = mean(vals); %#ok<AGROW>
            end
            level_labels{end+1} = sprintf('L%d', level); %#ok<AGROW>
        end

        bar(level_means, 'BarWidth', 0.6);
        xticks(1:numel(level_labels));
        xticklabels(level_labels);
        ylabel('Mean Density', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('%s System', sys), 'FontSize', 12, 'FontWeight', 'bold');
        ylim([0 1]);
        grid on;
        set(gca, 'GridAlpha', 0.3);
    end

    sgtitle('Neurotransmitter Systems Across CHARM Levels', 'FontSize', 16, 'FontWeight', 'bold');

    png_file = fullfile(output_dir, 'System_Comparison.png');
    pdf_file = fullfile(output_dir, 'System_Comparison.pdf');
    try
        exportgraphics(f, png_file, 'Resolution', 300);
        exportgraphics(f, pdf_file, 'ContentType', 'vector');
    catch
        saveas(f, png_file);
    end
    close(f);
end


function create_hierarchy_overview(all_levels_data, receptor_names, output_dir)
    f = figure('Position', [50 50 1400 800], 'Color', 'w');

    num_receptors = numel(receptor_names);
    level_receptor_matrix = NaN(5, num_receptors);
    row_labels = {};
    row_idx = 0;

    for level = 2:6
        if numel(all_levels_data) < level
            continue;
        end
        if isempty(all_levels_data(level).values)
            continue;
        end

        row_idx = row_idx + 1;
        for r = 1:num_receptors
            vals = all_levels_data(level).values(:, r);
            vals = vals(~isnan(vals));
            if isempty(vals)
                level_receptor_matrix(row_idx, r) = 0;
            else
                level_receptor_matrix(row_idx, r) = mean(vals);
            end
        end
        row_labels{row_idx} = sprintf('Level %d', level); %#ok<AGROW>
    end

    level_receptor_matrix(isnan(level_receptor_matrix)) = 0;

    display_receptor_names = { ...
        'AMPA (Glutamate)', 'Kainate (Glutamate)', 'NMDA (Glutamate)', ...
        'GABA_A', 'GABA_{ABZ}', 'GABA_B', ...
        'M_1 (ACh)', 'M_2 (ACh)', 'M_3 (ACh)', ...
        '5-HT_{1A} (Serotonin)', '5-HT_{2A} (Serotonin)', ...
        '\alpha_1 (Noradrenaline)', '\alpha_2 (Noradrenaline)', ...
        'D_1 (Dopamine)'};

    imagesc(level_receptor_matrix);
    colormap(parula);
    cb = colorbar;
    cb.Label.String = 'Mean Normalized Density';
    cb.FontSize = 12;

    xticks(1:num_receptors);
    xticklabels(display_receptor_names);
    xtickangle(45);

    yticks(1:numel(row_labels));
    yticklabels(row_labels);

    ylabel('CHARM Hierarchy Level', 'FontSize', 13, 'FontWeight', 'bold');
    xlabel('Receptor Type', 'FontSize', 13, 'FontWeight', 'bold');
    title('Receptor Density Across Hierarchical Levels (Regional Averages)', ...
        'FontSize', 14, 'FontWeight', 'bold');

    set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'tex');
    axis tight;

    png_file = fullfile(output_dir, 'Hierarchy_Overview.png');
    pdf_file = fullfile(output_dir, 'Hierarchy_Overview.pdf');
    try
        exportgraphics(f, png_file, 'Resolution', 300);
        exportgraphics(f, pdf_file, 'ContentType', 'vector');
    catch
        saveas(f, png_file);
    end
    close(f);
end
