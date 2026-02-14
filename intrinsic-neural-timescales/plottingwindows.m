figure;

for w_idx = 1:length(window_lengths)
    int_data = data_monks(w_idx:length(window_lengths):end);
    all_ints = cellfun(@(x) x.int, int_data, 'UniformOutput', false); % Extract INTs
    all_ints = cat(2, all_ints{:}); % Concatenate across runs

    % Compute average and error bars
    avg_int = mean(all_ints, 1, 'omitnan');
    std_int = std(all_ints, 0, 1, 'omitnan');
    sem_int = std_int / sqrt(length(monkeys)); % Standard Error of the Mean

    % Compute time points for each window center
    step_size = round(window_sizes(w_idx) / 2);
    num_windows = size(all_ints, 2);
    time_points = (step_size:step_size:(num_windows * step_size)) * TR; % Convert to seconds

    % Downsample for cleaner plot
    downsample_factor = 2; % Plot every 2nd point
    time_points = time_points(1:downsample_factor:end);
    avg_int = avg_int(1:downsample_factor:end);
    sem_int = sem_int(1:downsample_factor:end);

    % Create a new subplot for each window length
    subplot(length(window_lengths), 1, w_idx);
    hold on;
    errorbar(time_points, avg_int, sem_int, 'o-', ...
         'LineWidth', 1.5, 'Color', colors(w_idx), ...
         'MarkerSize', 3, 'CapSize', 8, ...
         'DisplayName', [num2str(window_lengths(w_idx)) ' seconds']);
    grid on; % Subtle grid for readability
    grid minor;
    xlabel('Time (seconds)', 'FontSize', 10);
    ylabel('Average INT', 'FontSize', 10);
    title(['Temporal Dynamics for ', num2str(window_lengths(w_idx)), 's Window'], 'FontSize', 12);
    legend('Location', 'best', 'FontSize', 8); % Compact legend
    hold off;

    % Adjust subplot spacing
    set(gca, 'Box', 'off', 'XMinorGrid', 'on', 'YMinorGrid', 'on'); % Add minor grid lines
end

% Improve figure layout
set(gcf, 'Color', 'w');
tight_layout = get(gcf, 'Position');
tight_layout(4) = tight_layout(4) * 1.1; % Increase vertical space
set(gcf, 'Position', tight_layout);
