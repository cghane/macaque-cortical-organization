% Iterates over monkeys and runs to compute timescales from processed fMRI data
% Initialize variables
data_monks = [];

% Parameters
no_runs = 5;
base_path = '/mnt/scratch/NHP4CYRUS/data_dump/';
folder_info = dir(base_path);
folder_names = {folder_info([folder_info.isdir]).name}; % Only directories
monkeys = folder_names(~ismember(folder_names, {'Aragorn', 'Andrea', '.', '..'})); % Exclude 'Andrea', '.' and '..'

n = 1; % Counter for storing data
for i = 1:length(monkeys)
    monkey = monkeys{i};
    for area = {'cortical'}
        area_name = area{1}; % Extract area name as string
        
        % Preallocate based on the maximum expected number of voxels
        num_voxels = 0;
        for run = 1:no_runs
            file_path = fullfile(base_path, monkey, area_name, ...
                                 ['ROIdump_', area_name, num2str(run), '.1D']);
            if isfile(file_path)
                data = importdata(file_path);
                num_voxels = max(num_voxels, size(data, 1));
            end
        end
        
        % Initialize results for this monkey and area
        int = zeros(num_voxels, no_runs); % Initialize int with appropriate size
        ac_magnitude = zeros(num_voxels, no_runs); % Initialize ac_magnitude
        lags = zeros(num_voxels, no_runs); % Initialize lags
        
        for run = 1:no_runs
            % Construct path dynamically
            file_path = fullfile(base_path, monkey, area_name, ...
                                 ['ROIdump_', area_name, num2str(run), '.1D']);
            if ~isfile(file_path)
                continue;
            end
            
            % Load data
            try
                data = importdata(file_path);
            catch
                continue;
            end
            
            % Initialize variables for this run
            num_voxels = size(data, 1);
            acf_data = zeros(num_voxels, size(data, 2));
            sum_acf = zeros(num_voxels, 1);
            lag_acf = zeros(num_voxels, 1);
            
            % Compute autocorrelation and timescales
            for ii = 1:num_voxels
                voxel = data(ii, :);
                
                % Calculate maximum possible lag length
                max_lag = length(voxel) - 1;
                
                try
                    % Compute autocorrelation with the correct lag
                    acf = autocorr(voxel, max_lag);
                    acf_data(ii, 1:max_lag + 1) = acf;
                    
                    % Compute INT as sum of positive ACF values
                    positive_acf = acf > 0;
                    sum_acf(ii) = sum(acf(positive_acf));
                    
                    % Find the first lag where ACF crosses zero
                    lag_acf(ii) = find(~positive_acf, 1) - 1;
                    if isempty(lag_acf(ii))
                        lag_acf(ii) = max_lag; % If no zero-crossing, use max lag
                    end
                catch err
                    fprintf('Error computing autocorr for voxel %d: %s\n', ii, err.message);
                    continue;
                end
            end
            
            % Store results for each run
            int(1:num_voxels, run) = sum_acf;
            ac_magnitude(1:num_voxels, run) = acf_data(:, 2); % Lag-one ACF value
            lags(1:num_voxels, run) = lag_acf;
        end
        
        % Convert lags to time (e.g., with 1.1 ms per lag step, if needed)
        lags_time = lags * 1.1; % Adjust the multiplier based on your time resolution
        
        % Store results for this monkey and area
        timescales.monkey = monkey;
        timescales.area = area_name;
        timescales.lags = lags;
        timescales.ac_mag = ac_magnitude;
        timescales.time = lags_time;
        timescales.int = int;
        
        data_monks{n} = timescales;
        n = n + 1;
    end
end
