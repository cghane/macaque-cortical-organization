% Loads each run’s ROI time series, chops it into 30/60/90s chunks, computes per-voxel intrinsic timescale (sum of positive ACF) for each chunk
% Then, averages across chunks/runs, stores results, then plots how the average timescale changes with window length.
% Modified Code to Compute INTs for 30, 60, and 90-Second Windows
int = []; % Intrinsic neural timescales
data_monks = [];

% Parameters
no_runs = 5;
n = 1;
base_path = '/mnt/scratch/NHP4CYRUS/data_dump/';
folder_info = dir(base_path);
folder_names = {folder_info([folder_info.isdir]).name}; % Only directories
%monkeys = folder_names(~ismember(folder_names, {'.', '..'})); % Exclude '.' and '..'
monkeys = {'Bilbo'}
% Define window lengths in seconds
window_lengths = [30, 60, 90];
TR = 1.1; % Repetition time in seconds
window_sizes = round(window_lengths / TR); % Convert to data points

for i = 1:length(monkeys)
    monkey = monkeys{i};
    for area = {'cortical'}
        area_name = area{1};
        int_by_window = struct();

        for w_idx = 1:length(window_sizes)
            window_size = window_sizes(w_idx);
            int_by_window(w_idx).window_length = window_lengths(w_idx);
            int_by_window(w_idx).int = [];
        end

        for run = 1:no_runs
            try
                data = importdata([base_path, monkey, '/', area_name, '/ROIdump_', area_name, num2str(run), '.1D']);
            catch
                fprintf('Skipping run %d for monkey %s: data load failed.\n', run, monkey);
                continue
            end

            if isempty(data) || size(data, 2) < max(window_sizes)
                fprintf('Skipping monkey %s, run %d: insufficient data.\n', monkey, run);
                continue;
            end

            % Segment data into windows and compute ACF for each segment
            for w_idx = 1:length(window_sizes)
                window_size = window_sizes(w_idx);
                num_windows = floor(size(data, 2) / window_size);
                int_window = zeros(size(data, 1), num_windows);

                for win = 1:num_windows
                    start_idx = (win - 1) * window_size + 1;
                    end_idx = start_idx + window_size - 1;
                    segment = data(:, start_idx:end_idx);

                    if size(segment, 2) < 2
                        fprintf('Skipping window %d: insufficient time points (%d).\n', win, size(segment, 2));
                        continue;
                    end

                    acf_data = zeros(size(segment, 1), window_size);
                    for voxel = 1:size(segment, 1)
                        % Check for zero variance
                        if var(segment(voxel, :)) == 0
                            fprintf('Skipping voxel %d in window %d: zero variance.\n', voxel, win);
                            continue;
                        end

                        % Dynamically set max_lag
                        max_lag = min(30, size(segment, 2) - 1);
                        if max_lag < 1
                            fprintf('Skipping voxel %d: insufficient lags (%d).\n', voxel, max_lag);
                            continue;
                        end

                        % Compute ACF
                        acf = autocorr(segment(voxel, :), max_lag);
                        acf_data(voxel, 1:max_lag+1) = acf;
                    end

                    % Compute INT for each voxel in this window
                    sum_acf = zeros(size(acf_data, 1), 1);
                    for voxel = 1:size(acf_data, 1)
                        for lag = 1:size(acf_data, 2)
                            if acf_data(voxel, lag) > 0
                                sum_acf(voxel) = sum_acf(voxel) + acf_data(voxel, lag);
                            else
                                break;
                            end
                        end
                    end

                    int_window(:, win) = sum_acf; % Store INTs for this window
                end

                % Average INTs across all windows for this length
                int_by_window(w_idx).int = [int_by_window(w_idx).int, mean(int_window, 2, 'omitnan')];
            end
        end

        % Store Results
        for w_idx = 1:length(window_sizes)
            timescales.monkey = monkey;
            timescales.area = area;
            timescales.window_length = int_by_window(w_idx).window_length;
            timescales.int = int_by_window(w_idx).int;
            data_monks{n} = timescales;
            n = n + 1;
        end
    end
end

% Plot Results w/ Error Bars
figure;
colors = ['r', 'g', 'b'];
hold on;

for w_idx = 1:length(window_lengths)
    % Extract data for the current window length
    int_data = data_monks(w_idx:length(window_lengths):end);
    avg_int = cellfun(@(x) mean(x.int, 'all', 'omitnan'), int_data);
    std_int = cellfun(@(x) std(x.int, 0, 'all', 'omitnan'), int_data);
    sem_int = std_int / sqrt(length(monkeys)); % Standard Error of the Mean

    % Plot with error bars
    errorbar(window_lengths(w_idx) * ones(size(avg_int)), avg_int, sem_int, ...
             'o-', 'Color', colors(w_idx), 'DisplayName', [num2str(window_lengths(w_idx)) ' seconds']);
end

xlabel('Window Length (seconds)');
ylabel('Average INT');
title('Intrinsic Neural Timescales Across Window Lengths');
legend('show');
hold off;
