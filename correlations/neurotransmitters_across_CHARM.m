%% computes the average value of each neurotransmitter across ROIs for CHARM levels C2–C6 and generates a saved bar chart comparing those averages across levels.
%Define file prefix
prefix = '/mnt/scratch/NHP4CYRUS/';
output_dir = [prefix 'output/'];

% Ensure output directory exists
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Neurotransmitter names
neuro_names = {'AMPA', 'Kainate', 'NMDA', 'GABA_A', 'GABA_A/BZ', 'GABA_B', ...
               'M_1', 'M_2', 'M_3', '5-HT_1A', '5-HT_2A', 'α_1', 'α_2', 'D_1'};

% Load RM_dump_nan.1D for neurotransmitter data
try
    neuro_data = load([prefix 'RM_dump_nan.1D']); % 17 columns: XYZ (1-3), neurotransmitters (4-17)
catch e
    error('Failed to load RM_dump_nan.1D: %s', e.message);
end
if size(neuro_data, 2) ~= 17
    error('RM_dump_nan.1D should have 17 columns, found %d.', size(neuro_data, 2));
end

% Initialize storage for average neurotransmitter values per charm level
charm_levels = 2:6;
num_neuro = length(neuro_names); % 14 neurotransmitters
avg_neuro_by_charm = zeros(length(charm_levels), num_neuro); % Rows: charm levels, Columns: neurotransmitters
valid_data = false(length(charm_levels), 1); % Track if data exists for each charm level

% Loop over each charm level (C2 to C6)
for charm_num = charm_levels
    charm_idx = charm_num - 1; % Map charm 2->1, 3->2, ..., 6->5 for indexing
    atlas_file = [prefix sprintf('/1dfiles/C%d_whole_brain.1D', charm_num)];
    
    % Load atlas file for the charm level
    try
        atlas_data = importdata(atlas_file); % Atlas labels
    catch e
        warning('Failed to load %s: %s. Skipping charm level %d.', atlas_file, e.message, charm_num);
        continue;
    end
    if size(atlas_data, 2) ~= 1
        error('%s should have 1 column, found %d.', atlas_file, size(atlas_data, 2));
    end
    atlas = atlas_data; % Atlas labels
    
    % Get all unique ROI labels in the atlas
    roi_labels = unique(atlas);
    roi_labels = roi_labels(roi_labels > 0); % Exclude background (label 0, if present)
    
    % Initialize storage for neurotransmitter values for this charm level
    neuro_values_per_roi = cell(1, num_neuro); % Store values for each neurotransmitter
    for i = 1:num_neuro
        neuro_values_per_roi{i} = [];
    end
    
    % Loop over all ROI labels in this atlas
    for label = roi_labels'
        % Find voxel indices for the current ROI
        crd_acc = find(atlas == label);
        if isempty(crd_acc)
            warning('No voxels found for label %d in %s. Skipping.', label, atlas_file);
            continue;
        end
        
        % Filter neuro_data for ROI voxels
        if max(crd_acc) > size(neuro_data, 1)
            warning('crd_acc indices for label %d in %s exceed rows in RM_dump_nan.1D (%d). Skipping.', ...
                    label, atlas_file, size(neuro_data, 1));
            continue;
        end
        data_acc = neuro_data(crd_acc, 4:end); % Neurotransmitter values (columns 4-17)
        
        % Check for NaN or Inf
        if any(isnan(data_acc(:))) || any(isinf(data_acc(:)))
            warning('NaN or Inf values in neurotransmitter data for label %d in %s. Skipping invalid values.', ...
                    label, atlas_file);
            data_acc = data_acc(~any(isnan(data_acc) | isinf(data_acc), 2), :);
            if isempty(data_acc)
                warning('No valid neurotransmitter values for label %d in %s after removing NaN/Inf. Skipping.', ...
                        label, atlas_file);
                continue;
            end
        end
        
        % Compute average for each neurotransmitter in this ROI
        roi_avg = mean(data_acc, 1); % Average across voxels for each of the 14 neurotransmitters
        
        % Append averages to the respective neurotransmitter's collection
        for i = 1:num_neuro
            neuro_values_per_roi{i} = [neuro_values_per_roi{i}; roi_avg(i)];
        end
    end
    
    % Compute average for each neurotransmitter across all ROIs in this charm level
    for i = 1:num_neuro
        if ~isempty(neuro_values_per_roi{i})
            avg_neuro_by_charm(charm_idx, i) = mean(neuro_values_per_roi{i});
            valid_data(charm_idx) = true;
        else
            avg_neuro_by_charm(charm_idx, i) = NaN; % No data for this neurotransmitter
        end
    end
end

% Check for charm levels with no data
for i = 1:length(charm_levels)
    if ~valid_data(i)
        warning('No valid neurotransmitter data available for charm level %d.', charm_levels(i));
    end
end

% Create bar chart
figure('Position', [100, 100, 1200, 800]);
bar(1:num_neuro, avg_neuro_by_charm');
set(gca, 'XTick', 1:num_neuro, 'XTickLabel', neuro_names, 'XTickLabelRotation', 45);
xlabel('Neurotransmitter');
ylabel('Average Value');
title('Average Neurotransmitter Values Across All ROIs by Charm Level');
legend(arrayfun(@(x) sprintf('Charm Level %d', x), charm_levels, 'UniformOutput', false), ...
       'Location', 'best');
grid on;
set(gcf, 'Color', 'white');

% Save the chart
output_png = [output_dir 'avg_neuro_per_transmitter_by_charm.png'];
try
    saveas(gcf, output_png);
    fprintf('Saved bar chart to %s\n', output_png);
catch e
    warning('Failed to save bar chart: %s', e.message);
end

% Close figure to save memory
close(gcf);
