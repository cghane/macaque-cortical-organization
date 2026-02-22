%% computes voxel-wise correlations between INT (four segments) and 14 neurotransmitters within each ROI, averages those correlations across regions, saves the results, and visualizes them with heatmaps and scatter plots.
% Step 1: Load and Extract Neurotransmitter Data
% Load the neurotransmitter data
RM_dump = importdata('/mnt/scratch/NHP4CYRUS/RM_dump.1D');

% Debug: Check RM_dump size
fprintf('Size of RM_dump: %s\n', mat2str(size(RM_dump)));

% Load the XYZ coordinates file
ijk = importdata('/mnt/scratch/NHP4CYRUS/ijk_whole_brain.1D'); % Contains XYZ coordinates for all atlas voxels

% Debug: Check size and sample of ijk
fprintf('Size of ijk_whole_brain.1D: %s\n', mat2str(size(ijk)));
fprintf('Sample ijk (first 5 rows):\n');
disp(ijk(1:5, :));
fprintf('ijk XYZ range (min, max):\n');
disp([min(ijk(:, 1:3)); max(ijk(:, 1:3))]);

% Define ROIs, their atlas labels, and corresponding atlas files
rois = struct();
rois.dlPFC = struct('label', 54, 'atlas', 'C3');
rois.latOFC = struct('label', 25, 'atlas', 'C3');
rois.PMd = struct('label', 81, 'atlas', 'C5');
rois.ACC = struct('label', 3, 'atlas', 'C3');
rois.OFC = struct('label', 16, 'atlas', 'C2');
rois.LPFC = struct('label', 50, 'atlas', 'C2');
rois.LIP = struct('label', 113, 'atlas', 'C4');
rois.MT = struct('label', 232, 'atlas', 'C5');
rois.S2 = struct('label', 95, 'atlas', 'C4');
rois.S1 = struct('label', 92, 'atlas', 'C4');

% List of ROI names for iteration
roi_names = fieldnames(rois);

% Neurotransmitter names (for columns 4–17)
neurotransmitters = {'AMPA', 'Kainate', 'NMDA', 'GABA_A', 'GABA_A_BZ', 'GABA_B', ...
                     'M_1', 'M_2', 'M_3', '5-HT_1A', '5-HT_2A', 'alpha_1', 'alpha_2', 'D_1'};

% Structure to store neurotransmitter data for each ROI
neurotransmitter_data = struct();
% Structure to store INT data for each ROI
int_data = struct();

% Extract data for each ROI using the correct atlas
for i = 1:length(roi_names)
    roi = roi_names{i};
    label = rois.(roi).label;
    atlas_file = ['/mnt/scratch/NHP4CYRUS/1dfiles/' rois.(roi).atlas, '_whole_brain.1D']; % Add full path
    
    % Load the appropriate atlas
    atlas = importdata(atlas_file);
    
    % Find coordinates in the atlas for this ROI
    crd_acc = find(atlas == label);
    
    % Skip if no coordinates found
    if isempty(crd_acc)
        fprintf('Error: Region %s has no coordinates in atlas %s. Skipping.\n', ...
                roi, atlas_file);
        continue;
    end
    
    % Validate crd_acc indices for ijk
    if max(crd_acc) > size(ijk, 1)
        fprintf('Error: Region %s has invalid coordinates for ijk (max crd: %d, ijk rows: %d). Skipping.\n', ...
                roi, max(crd_acc), size(ijk, 1));
        continue;
    end
    
    % Get XYZ coordinates for the ROI from ijk (first 3 columns)
    roi_xyz = ijk(crd_acc, 1:3); % Ensure 3 columns (X, Y, Z)
    fprintf('Region: %s, roi_xyz size: %s\n', roi, mat2str(size(roi_xyz)));
    fprintf('Region: %s, roi_xyz range (min, max):\n', roi);
    disp([min(roi_xyz); max(roi_xyz)]);
    
    % Load one INT file to get the XYZ coordinates and INT data
    file_name = sprintf('%s_avg_segment1.mat', roi);
    try
        data = importdata(file_name);
        if isstruct(data) && isfield(data, 'data_acc')
            int_data_matrix = data.data_acc; % 13581 x 4 matrix
        else
            int_data_matrix = data; % 13581 x 4 matrix
        end
    catch
        fprintf('Error: Could not load file %s. Skipping region.\n', file_name);
        continue;
    end
    
    % Validate INT data matrix
    if size(int_data_matrix, 2) < 4
        fprintf('Error: Region %s INT data has %d columns, expected at least 4. Skipping.\n', ...
                roi, size(int_data_matrix, 2));
        continue;
    end
    fprintf('Region: %s, int_data_matrix size: %s\n', roi, mat2str(size(int_data_matrix)));
    
    % Extract XYZ coordinates from INT data (columns 1–3)
    int_xyz = int_data_matrix(:, 1:3); % 13581 x 3 matrix
    fprintf('Region: %s, int_xyz size: %s\n', roi, mat2str(size(int_xyz)));
    fprintf('Region: %s, int_xyz range (min, max):\n', roi);
    disp([min(int_xyz); max(int_xyz)]);

    % Debug: Sample XYZ coordinates
    fprintf('Region: %s, Sample roi_xyz (first 5 rows):\n', roi);
    disp(roi_xyz(1:min(5, size(roi_xyz, 1)), :));
    fprintf('Region: %s, Sample int_xyz (first 5 rows):\n', roi);
    disp(int_xyz(1:5, :));
    
    % Compute dynamic shifts based on median differences to handle outliers
    roi_xyz_median = median(roi_xyz, 1);
    int_xyz_median = median(int_xyz, 1);
    shifts = roi_xyz_median - int_xyz_median; % [X_shift, Y_shift, Z_shift]
    fprintf('Region: %s, Computed shifts (X, Y, Z): %.2f, %.2f, %.2f\n', ...
            roi, shifts(1), shifts(2), shifts(3));
    
    % Apply shifts to roi_xyz
    roi_xyz_shifted = roi_xyz - repmat(shifts, size(roi_xyz, 1), 1); % Use repmat for compatibility
    
    % Find matching rows in INT data
    % Version-compatible approach for ismembertol
    tolerance = 1; % Adjust as needed
    try
        % Try with 'ByRows' parameter (newer MATLAB)
        [is_member, matched_indices] = ismembertol(int_xyz, roi_xyz_shifted, tolerance, 'ByRows', true);
    catch
        % Fallback for older MATLAB versions
        is_member = false(size(int_xyz, 1), 1);
        matched_indices = zeros(size(int_xyz, 1), 1);
        
        % Manual row-by-row comparison within tolerance
        for row = 1:size(int_xyz, 1)
            diffs = roi_xyz_shifted - repmat(int_xyz(row,:), size(roi_xyz_shifted, 1), 1);
            maxDiffs = max(abs(diffs), [], 2);
            minDiff = min(maxDiffs);
            if minDiff <= tolerance
                is_member(row) = true;
                [~, idx] = min(maxDiffs);
                matched_indices(row) = idx;
            end
        end
    end
    
    % Fallback: Nearest-neighbor matching if no exact matches
    if sum(is_member) < 0.1 * length(crd_acc) % Less than 10% of expected matches
        fprintf('Warning: Region %s has few matches (%d). Trying nearest-neighbor matching.\n', ...
                roi, sum(is_member));
        
        % Version-compatible approach for pdist2
        try
            % Try with pdist2 (newer MATLAB)
            distances = pdist2(int_xyz, roi_xyz_shifted, 'euclidean');
        catch
            % Fallback for older MATLAB versions - manual distance calculation
            distances = zeros(size(int_xyz, 1), size(roi_xyz_shifted, 1));
            for ii = 1:size(int_xyz, 1)
                for jj = 1:size(roi_xyz_shifted, 1)
                    distances(ii,jj) = sqrt(sum((int_xyz(ii,:) - roi_xyz_shifted(jj,:)).^2));
                end
            end
        end
        
        [~, matched_indices] = min(distances, [], 2); % Closest roi_xyz for each int_xyz
        valid_indices = 1:size(int_xyz, 1); % Use all int_xyz indices
        matched_roi_indices = matched_indices; % Corresponding roi_xyz indices
    else
        valid_indices = find(is_member); % Indices in int_xyz that match
        matched_roi_indices = matched_indices(is_member); % Corresponding roi_xyz indices
    end
    
    % Warn if no matches or fewer matches than expected
    if isempty(valid_indices)
        fprintf('Error: Region %s has no matching XYZ coordinates in INT data (tolerance: %.2f). Skipping.\n', ...
                roi, tolerance);
        continue;
    end
    if length(valid_indices) < length(crd_acc)
        fprintf('Warning: Region %s had %d atlas coordinates, but only %d matched INT data XYZ coordinates (tolerance: %.2f).\n', ...
                roi, length(crd_acc), length(valid_indices), tolerance);
    end
    
    % Validate indices for RM_dump
    if max(valid_indices) > size(RM_dump, 1)
        fprintf('Error: Region %s has invalid indices for RM_dump (max idx: %d, RM_dump rows: %d). Skipping.\n', ...
                roi, max(valid_indices), size(RM_dump, 1));
        continue;
    end
    
    % Extract neurotransmitter data for these indices
    neurotransmitter_data.(roi) = RM_dump(valid_indices, 4:end); % Voxel-level data
    fprintf('Region: %s, Neurotransmitter voxels: %d\n', roi, length(valid_indices));
    
    % Load INT data for each segment
    for seg = 1:4
        % Construct the file name (e.g., 'dlPFC_avg_segment1.mat')
        file_name = sprintf('%s_avg_segment%d.mat', roi, seg);
        
        % Load the INT data
        try
            data = importdata(file_name);
        catch
            fprintf('Error: Could not load file %s. Skipping segment.\n', file_name);
            continue;
        end
        
        % Extract INT values from column 4 for matched indices
        if isstruct(data) && isfield(data, 'data_acc')
            if size(data.data_acc, 2) < 4
                fprintf('Error: Region %s, Segment %d INT data has %d columns, expected at least 4. Skipping segment.\n', ...
                        roi, seg, size(data.data_acc, 2));
                continue;
            end
            int_values_segment = data.data_acc(valid_indices, 4);
        else
            if size(data, 2) < 4
                fprintf('Error: Region %s, Segment %d INT data has %d columns, expected at least 4. Skipping segment.\n', ...
                        roi, seg, size(data, 2));
                continue;
            end
            int_values_segment = data(valid_indices, 4);
        end
        
        % Store INT data
        int_data.(roi).(['segment', num2str(seg)]) = int_values_segment;
        fprintf('Region: %s, Segment %d, INT voxels: %d, any NaN: %d, variance: %.4f\n', ...
                roi, seg, length(int_values_segment), any(isnan(int_values_segment)), var(int_values_segment));
    end
end

% Step 2: Correlate INT with Neurotransmitter Data (Voxel-wise)
% Number of ROIs
n_rois = length(roi_names);

% Number of neurotransmitters
n_neurotransmitters = length(neurotransmitters);

% Matrix to store correlations: 4 segments x 14 neurotransmitters
correlations = zeros(4, n_neurotransmitters);

% Compute correlations for each segment and neurotransmitter
for seg = 1:4
    for nt = 1:n_neurotransmitters
        % Collect correlations across all regions
        region_corrs = zeros(n_rois, 1);
        for i = 1:n_rois
            roi = roi_names{i};
            
            % Check if data exists for this region
            if ~isfield(int_data, roi) || ~isfield(int_data.(roi), ['segment', num2str(seg)]) || isempty(int_data.(roi).(['segment', num2str(seg)]))
                region_corrs(i) = NaN;
                continue;
            end
            
            % INT values for this segment and region (voxel-level)
            seg_field = ['segment', num2str(seg)];
            int_vals = int_data.(roi).(seg_field);
            
            % Neurotransmitter values for this neurotransmitter and region (voxel-level)
            nt_vals = neurotransmitter_data.(roi)(:, nt);
            
            % Ensure same dimensions
            if length(int_vals) ~= length(nt_vals)
                fprintf('Warning: Dimension mismatch for %s, segment %d, NT %s. Skipping.\n', ...
                        roi, seg, neurotransmitters{nt});
                region_corrs(i) = NaN;
                continue;
            end
            
            % Debug output
            fprintf('Region: %s, Segment: %d, Neurotransmitter: %s\n', roi, seg, neurotransmitters{nt});
            fprintf('INT size: %d, NT size: %d\n', length(int_vals), length(nt_vals));
            fprintf('INT variance: %.4f, NT variance: %.4f\n', var(int_vals), var(nt_vals));
            fprintf('INT NaN: %d, NT NaN: %d\n', sum(isnan(int_vals)), sum(isnan(nt_vals)));
            
            % Remove NaN values before correlation
            valid_idx = ~isnan(int_vals) & ~isnan(nt_vals);
            int_vals_clean = int_vals(valid_idx);
            nt_vals_clean = nt_vals(valid_idx);
            
            % Compute voxel-wise correlation
            if length(int_vals_clean) > 1 && var(int_vals_clean) > 0 && var(nt_vals_clean) > 0
                r = corrcoef(int_vals_clean, nt_vals_clean);
                if size(r, 1) >= 2 && size(r, 2) >= 2
                    region_corrs(i) = r(1, 2);
                else
                    region_corrs(i) = NaN;
                    fprintf('Correlation resulted in invalid matrix size: %s\n', mat2str(size(r)));
                end
            else
                region_corrs(i) = NaN;
                fprintf('Correlation skipped due to invalid data (too few voxels, zero variance, or NaNs)\n');
            end
        end
        % Average correlation across regions
        correlations(seg, nt) = nanmean(region_corrs);
        fprintf('Correlation for Segment %d, Neurotransmitter %s: %.4f\n', seg, neurotransmitters{nt}, correlations(seg, nt));
    end
end

% Step 2.5: Save Results to a .mat File
% Save correlations, neurotransmitter_data, and int_data to a descriptive .mat file
mat_filename = 'INT_Neurotransmitter_Correlations_2025.mat';
save(mat_filename, 'correlations', 'neurotransmitter_data', 'int_data', 'neurotransmitters', 'roi_names');
fprintf('Results saved to %s\n', mat_filename);

% Step 3: Visualize the Results
% Check for version compatibility with various functions
matlab_version = version;
fprintf('MATLAB Version: %s\n', matlab_version);

% 3.1 Heatmap of Correlations
figure('Name', 'Correlation Heatmap', 'Position', [100, 100, 1000, 600]);

% Try the newer heatmap function first
try
    h = heatmap(neurotransmitters, {'Segment 1', 'Segment 2', 'Segment 3', 'Segment 4'}, correlations, ...
                'ColorbarVisible', true, ...
                'CellLabelColor', 'black', ...
                'CellLabelFormat', '%.2f', ...
                'FontSize', 12);
    
    % Set title and labels directly for version compatibility
    title('Correlation between INT Segments and Neurotransmitters', 'FontSize', 18);
    xlabel('Neurotransmitters', 'FontSize', 14);
    ylabel('INT Segments', 'FontSize', 14);
    
    % Try to set colormap - use parula if available, otherwise jet
    try
        colormap(parula);
    catch
        % Fallback to jet colormap for older MATLAB versions
        colormap(jet);
    end
catch
    % Fallback for very old MATLAB versions without heatmap function
    imagesc(correlations);
    colorbar;
    title('Correlation between INT Segments and Neurotransmitters', 'FontSize', 18);
    xlabel('Neurotransmitters', 'FontSize', 14);
    ylabel('INT Segments', 'FontSize', 14);
    set(gca, 'XTick', 1:length(neurotransmitters), 'XTickLabel', neurotransmitters, ...
        'YTick', 1:4, 'YTickLabel', {'Segment 1', 'Segment 2', 'Segment 3', 'Segment 4'}, ...
        'XTickLabelRotation', 45, 'FontSize', 12);
    
    % Try to set colormap - use parula if available, otherwise jet
    try
        colormap(parula);
    catch
        % Fallback to jet colormap for older MATLAB versions
        colormap(jet);
    end
end

% 3.2 Bar Plots for Each Segment
figure('Name', 'Correlation Bar Plots', 'Position', [100, 100, 1200, 900]);
for seg = 1:4
    subplot(2, 2, seg);
    bar(correlations(seg, :), 'FaceColor', [0.2, 0.6, 0.8]);
    set(gca, 'XTick', 1:n_neurotransmitters, 'XTickLabel', neurotransmitters, ...
            'XTickLabelRotation', 45, 'FontSize', 10);
    title(['Segment ', num2str(seg)], 'FontSize', 12);
    ylabel('Correlation Coefficient', 'FontSize', 11);
    ylim([-1, 1]);
    grid on;
    box on;
end

% Try using sgtitle (newer) or manual annotation (older)
try
    sgtitle('Correlation between INT Segments and Neurotransmitters', 'FontSize', 16);
catch
    % Fallback for older MATLAB versions
    annotation('textbox', [0.5, 0.95, 0, 0], 'String', 'Correlation between INT Segments and Neurotransmitters', ...
               'HorizontalAlignment', 'center', 'FontSize', 16, 'EdgeColor', 'none');
end

% 3.3 Scatter Plots for Significant Correlations
significant_threshold = 0.5;
[seg_idx, nt_idx] = find(abs(correlations) > significant_threshold);

for i = 1:length(seg_idx)
    seg = seg_idx(i);
    nt = nt_idx(i);
    
    % Collect data across regions for scatter plot
    all_int_vals = [];
    all_nt_vals = [];
    for j = 1:n_rois
        roi = roi_names{j};
        seg_field = ['segment', num2str(seg)];
        if isfield(int_data, roi) && isfield(int_data.(roi), seg_field)
            int_vals = int_data.(roi).(seg_field);
            nt_vals = neurotransmitter_data.(roi)(:, nt);
            
            % Ensure same dimensions and no NaNs
            if length(int_vals) == length(nt_vals)
                valid_idx = ~isnan(int_vals) & ~isnan(nt_vals);
                all_int_vals = [all_int_vals; int_vals(valid_idx)];
                all_nt_vals = [all_nt_vals; nt_vals(valid_idx)];
            end
        end
    end
    
    % Only create scatter plot if we have enough data points
    if length(all_int_vals) > 5
        figure('Name', ['Scatter Plot: Segment ', num2str(seg), ' vs ', neurotransmitters{nt}], ...
               'Position', [100, 100, 700, 500]);
        scatter(all_int_vals, all_nt_vals, 50, 'filled', 'MarkerFaceColor', [0.2, 0.6, 0.8], ...
                'MarkerFaceAlpha', 0.6);
        xlabel(['INT Segment ', num2str(seg)], 'FontSize', 12);
        ylabel(neurotransmitters{nt}, 'FontSize', 12);
        title(['INT Segment ', num2str(seg), ' vs ', neurotransmitters{nt}, ...
               ' (r = ', num2str(correlations(seg, nt), '%.2f'), ')'], 'FontSize', 14);
        grid on;
        set(gca, 'FontSize', 10);
    else
        fprintf('Not enough data points for scatter plot: Segment %d vs %s\n', seg, neurotransmitters{nt});
    end
end
