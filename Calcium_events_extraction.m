% Calcium Imaging Data Processing Script
% Purpose: This script processes calcium imaging data from a .mat file from Suite2P, applies
% filtering, detects spikes, bins the data, and calculates statistics such as
% frequency, amplitude, and inter-event intervals. Results are saved as .mat and
% .csv files for further analysis.
%
% Usage: Update the 'data_dir' variable with the path to your data directory.
% Ensure the 'Fall.mat' file from Suite2P is located in the specified directory structure.
%
% Dependencies: MATLAB Signal Processing Toolbox
%
% Author: [l.bijoch@nencki.edu.pl]
% Date: [02.07.2025]

function process_calcium_data()

    %% CONFIGURATION
    data_dir = 'C:\desktop'; % <- Update path
    mat_file = fullfile(data_dir, 'avg', 'suite2p', 'plane0', 'Fall.mat');
    bin_size = 15;  % Number of bins
    num_bins = 420; % Pre-defined number of bins
    frame_rate = 6300 / 1260; % Hz of recorded data
    example_cell = 14; % Number of trace to extract
    spike_threshold_factor = 0.4; % Define threshold method for spiking detection

    %% LOAD DATA
    data = load(mat_file);
    F = double(data.F(data.iscell(:, 1) == 1, :));

    %% DFF CALCULATION
    dF = (F - mean(F, 2)) ./ mean(F, 2);

    %% FILTERING
    dF = bandpass_filter(dF, 0.006, 0.55);

    %% SPIKE DETECTION
    [spike_indices, thresholded_dF] = detect_spikes(dF, spike_threshold_factor);

    %% VISUALIZATION
    visualize_cell_signal(dF, spike_indices, example_cell);

    %% SAVE THRESHOLDED DATA
    save(fullfile(data_dir, 'dFc.mat'), 'thresholded_dF');
    writematrix(thresholded_dF, fullfile(data_dir, 'dFc_filtered.csv'));

    %% BINNING
    [event_counts, event_amplitudes] = bin_data(thresholded_dF, bin_size, num_bins);
    save(fullfile(data_dir, 'event_number.mat'), 'event_counts');
    save(fullfile(data_dir, 'F_signal.mat'), 'event_amplitudes');

    %% STATISTICS
    [freq, amp, iei_stats] = compute_statistics(thresholded_dF, frame_rate);
    writematrix(freq, fullfile(data_dir, 'frequency.mat'));
    writematrix(amp, fullfile(data_dir, 'amplitude.mat'));
    writematrix(iei_stats.mean_intervals, fullfile(data_dir, 'interevent_interval.mat'));

    %% SAVE STRUCT
    processed_data = struct('dF', dF, 'thresholded_dF', thresholded_dF);
    save(fullfile(data_dir, 'processedData.mat'), 'processed_data');

    disp('Processing complete. Results saved.');

end

%% ---------- SUBFUNCTIONS ----------

function filtered = bandpass_filter(signal, high_cutoff, low_cutoff)
    [b_low, a_low] = butter(4, low_cutoff);
    [b_high, a_high] = butter(4, high_cutoff, 'high');
    filtered = zeros(size(signal));
    for i = 1:size(signal,1)
        tmp = filtfilt(b_low, a_low, signal(i,:));
        filtered(i,:) = filtfilt(b_high, a_high, tmp);
    end
end

function [spike_idx, thresholded] = detect_spikes(signal, threshold_factor)
    max_val = max(signal, [], 2);
    thresholds = max_val * threshold_factor;
    spike_idx = signal > thresholds;
    thresholded = signal;
    thresholded(~spike_idx) = 0;
end

function visualize_cell_signal(dF, spikes, cell_idx)
    figure;
    plot(dF(cell_idx,:), 'b'); hold on;
    scatter(find(spikes(cell_idx,:)), dF(cell_idx, spikes(cell_idx,:)), 'r', 'filled');
    title(['Cell ', num2str(cell_idx), ' dF/F and Spikes']);
    xlabel('Frames'); ylabel('dF/F');
    legend('Signal', 'Spikes'); hold off;
end

function [counts, amplitudes] = bin_data(signal, bin_size, num_bins)
    num_cells = size(signal, 1);
    counts = zeros(num_cells, num_bins);
    amplitudes = zeros(num_cells, num_bins);
    for bin = 1:num_bins
        idx = (bin-1)*bin_size + 1 : bin*bin_size;
        if max(idx) > size(signal,2), break; end
        bin_data = signal(:, idx);
        counts(:, bin) = sum(bin_data > 0, 2);
        amplitudes(:, bin) = sum(bin_data, 2);
    end
end

function [frequency, mean_amplitude, iei_stats] = compute_statistics(thresholded_dF, frame_rate)
    mask = thresholded_dF > 0;
    mean_amplitude = sum(thresholded_dF .* mask, 2) ./ sum(mask, 2);

    % Frequency per cell
    frequency = sum(mask, 2) / (size(thresholded_dF,2) / frame_rate);

    % Inter-event intervals
    num_cells = size(thresholded_dF,1);
    intervals = cell(num_cells,1);
    for i = 1:num_cells
        spikes = find(mask(i,:));
        intervals{i} = diff(spikes);
    end

    mean_intervals = cellfun(@(x) mean(x), intervals, 'UniformOutput', true);
    mean_interval_global = mean(mean_intervals(mean_intervals > 0));

    % Combined histogram for visualization (if needed)
    all_iei = vertcat(intervals{:});
    all_iei = all_iei(all_iei > 0);
    all_iei = sort(all_iei);
    num_bins = 20;
    bin_edges = round(linspace(1, numel(all_iei), num_bins+1));
    binned_means = arrayfun(@(i) mean(all_iei(bin_edges(i):bin_edges(i+1)-1)), 1:num_bins);

    iei_stats = struct( ...
        'mean_intervals', mean_intervals, ...
        'mean_interval_global', mean_interval_global, ...
        'binned_means', binned_means ...
    );
end