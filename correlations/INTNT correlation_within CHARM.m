%% loads voxelwise neurotransmitter maps; uses each CHARM atlas level (C2–C6) to compute ROI-wise mean neurotransmitter values
prefix = '/mnt/scratch/NHP4CYRUS/';
output_dir = [prefix 'output_per_roi_named/'];

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Neurotransmitter names
neuro_names = {'AMPA', 'Kainate', 'NMDA', 'GABA_A', 'GABA_A/BZ', 'GABA_B', ...
               'M_1', 'M_2', 'M_3', '5-HT_1A', '5-HT_2A', 'α_1', 'α_2', 'D_1'};

% Load neurotransmitter data
neuro_data = load([prefix 'RM_dump_nan.1D']);
if size(neuro_data, 2) ~= 17
    error('Expected 17 columns in RM_dump_nan.1D.');
end

% Load CHARM key table
roi_key_file = '/mnt/scratch/NHP4CYRUS/CHARM_key_table.csv';
roi_table = readtable(roi_key_file, 'Delimiter', '\t', 'ReadVariableNames', true);

% Create mappings for each charm level
num_levels = width(roi_table); % should match max charm level
level_maps = cell(1, num_levels);

for lvl = 1:num_levels
    raw_entries = roi_table{:, lvl};
    map_lvl = containers.Map('KeyType', 'double', 'ValueType', 'char');
    for i = 1:length(raw_entries)
        entry = raw_entries{i};
        if ~ischar(entry) && ~isstring(entry)
            continue;
        end
        parts = split(entry, ':');
        if numel(parts) ~= 2
            continue;
        end
        roi_num = str2double(strtrim(parts{1}));
        roi_name = strtrim(parts{2});
        if ~isnan(roi_num)
            map_lvl(roi_num) = roi_name;
        end
    end
    level_maps{lvl} = map_lvl;
end

% Process charm levels 1–6
for charm_num = 2:6
    atlas_file = [prefix sprintf('/1dfiles/C%d_whole_brain.1D', charm_num)];
    try
        atlas_data = importdata(atlas_file);
    catch e
        warning('Could not load atlas file for charm %d: %s', charm_num, e.message);
        continue;
    end
    if size(atlas_data, 2) ~= 1
        error('Atlas file %s should have 1 column.', atlas_file);
    end
    atlas = atlas_data;

    % Get ROI labels
    roi_labels = unique(atlas);
    roi_labels = roi_labels(roi_labels > 0);

    % Get mapping for this charm level
    roi_map = level_maps{charm_num};

    % Storage
    roi_names = cell(length(roi_labels), 1);
    roi_neuro_data = nan(length(roi_labels), length(neuro_names));

    for r = 1:length(roi_labels)
        label = roi_labels(r);
        if isKey(roi_map, label)
            roi_names{r} = roi_map(label);
        else
            roi_names{r} = sprintf('ROI_%d', label);
        end

        voxel_idx = find(atlas == label);
        roi_vals = neuro_data(voxel_idx, 4:end);
        roi_vals = roi_vals(~any(isnan(roi_vals) | isinf(roi_vals), 2), :);
        if isempty(roi_vals)
            continue;
        end
        roi_neuro_data(r, :) = mean(roi_vals, 1);
    end

    % Save CSV
    csv_file = [output_dir sprintf('charm%d_roi_neuro_named.csv', charm_num)];
    T = array2table(roi_neuro_data, 'VariableNames', neuro_names, 'RowNames', roi_names);
    writetable(T, csv_file, 'WriteRowNames', true);

    % Plot
    figure('Position', [100, 100, 1400, 800]);
    bar(roi_neuro_data, 'grouped');
    set(gca, 'XTick', 1:length(roi_labels), 'XTickLabel', roi_names, 'XTickLabelRotation', 45);
    xlabel('Region of Interest');
    ylabel('Average Neurotransmitter Value');
    title(sprintf('ROI-Level Neurotransmitter Values - Charm Level %d', charm_num));
    legend(neuro_names, 'Location', 'northeastoutside');
    grid on;
    set(gcf, 'Color', 'white');

    plot_file = [output_dir sprintf('charm%d_roi_neuro_named.png', charm_num)];
    saveas(gcf, plot_file);
    close(gcf);
end
