% Split RM_dump_nan.1D (aggregate file with neurotransmitter receptor distribution data) into 14 individual neurotransmitter files

% Define receptor names
receptor_names = {
    'AMPA_Glutamate', 'Kainate_Glutamate', 'NMDA_Glutamate', ...
    'GABA_A', 'GABA_A_BZ', 'GABA_B', ...
    'M1_Acetylcholine', 'M2_Acetylcholine', 'M3_Acetylcholine', ...
    'HT1A_Serotonin', 'HT2A_Serotonin', ...
    'Alpha1_Noradrenaline', 'Alpha2_Noradrenaline', ...
    'D1_Dopamine'
};

% Create output directory
output_dir = 'RM_dump_transmitters';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created directory: %s\n', output_dir);
end

% Load the data
fprintf('Loading RM_dump_nan.1D...\n');
data = readmatrix('RM_dump_nan.1D', 'FileType', 'text');

% Check data dimensions
fprintf('Data size: %d rows x %d columns\n', size(data, 1), size(data, 2));
if size(data, 2) < 17
    error('Expected at least 17 columns (3 coords + 14 transmitters), got %d', size(data, 2));
end

% Extract coordinates and transmitter values
coords = data(:, 1:3);            % xyz coordinates
transmitter_values = data(:, 4:17);  % 14 transmitter columns

% Split into individual transmitter files
for t = 1:14
    output_file = fullfile(output_dir, sprintf('%s.1D', receptor_names{t}));
    
    % Open file for writing
    fid = fopen(output_file, 'w');
    if fid == -1
        error('Could not create file: %s', output_file);
    end
    
    % Write header
    fprintf(fid, '# %s receptor data\n', receptor_names{t});
    fprintf(fid, '# Split from RM_dump_nan.1D\n');
    fprintf(fid, '# Format: x y z value\n');
    
    % Write data (x y z value format) - keep ALL rows including NaN
    valid_count = 0;
    nan_count = 0;
    for i = 1:size(coords, 1)
        value = transmitter_values(i, t);
        if isnan(value)
            fprintf(fid, '%.6f %.6f %.6f NaN\n', ...
                coords(i, 1), coords(i, 2), coords(i, 3));
            nan_count = nan_count + 1;
        else
            fprintf(fid, '%.6f %.6f %.6f %.6f\n', ...
                coords(i, 1), coords(i, 2), coords(i, 3), value);
            valid_count = valid_count + 1;
        end
    end
    fclose(fid);
    
    % Report progress
    total_rows = size(coords, 1);
    fprintf('Created %s: %d total rows (%d valid, %d NaN)\n', ...
        sprintf('%s.1D', receptor_names{t}), total_rows, valid_count, nan_count);
end

fprintf('\nSplit complete! Files saved in: %s\n', output_dir);
fprintf('Created %d transmitter files\n', length(receptor_names));

% List created files
fprintf('\nCreated files:\n');
files = dir(fullfile(output_dir, '*.1D'));
for i = 1:length(files)
    fprintf('  %s\n', files(i).name);
end
