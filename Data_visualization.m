%% Calcium Imaging Data Analysis and Visualization
% Author: l.bijoch@nencki.edu.pl | Date: 02.07.2025
% Description: Combines, normalizes, clusters, and visualizes DRD1 and DRD2 calcium imaging data.

clear;

%% 1. Parameters
blur_sigma = 2; % Define blur effect for generated heatmaps
colormap_name = 'seismic'; % fallback to MATLAB colormaps if slanCM not found
num_clusters = 3; % Define number of cluster e.g. based on dendrograms

% Define your directory with .mat files with extracted matrices with binned DRD1
% and DRD2 calcium transients
base_dir = fullfile('C:', 'spikes');
drd1_dir = fullfile(base_dir, 'DRD1');
drd2_dir = fullfile(base_dir, 'DRD2');
% Define your output directory
drd1_output_dir = fullfile(drd1_dir, 'output2');
drd2_output_dir = fullfile(drd2_dir, 'output2');
if ~exist(drd1_output_dir, 'dir'), mkdir(drd1_output_dir); end
if ~exist(drd2_output_dir, 'dir'), mkdir(drd2_output_dir); end

%% 2. Load Data
drd1_data = []; drd2_data = [];
drd1_mean_per_animal = []; drd2_mean_per_animal = [];

drd2_files = dir(fullfile(drd2_dir, '*.mat'));
assert(~isempty(drd2_files), 'No .mat files in DRD2 dir.');
for file = drd2_files'
    vars = load(fullfile(file.folder, file.name));
    for f = fieldnames(vars)'
        d = vars.(f{1});
        drd2_data = [drd2_data; d];
        drd2_mean_per_animal = [drd2_mean_per_animal; mean(d, 1)];
    end
end

drd1_files = dir(fullfile(drd1_dir, '*.mat'));
assert(~isempty(drd1_files), 'No .mat files in DRD1 dir.');
for file = drd1_files'
    vars = load(fullfile(file.folder, file.name));
    for f = fieldnames(vars)'
        d = vars.(f{1});
        drd1_data = [drd1_data; d];
        drd1_mean_per_animal = [drd1_mean_per_animal; mean(d, 1)];
    end
end

%% 3. Save Mean Activity
activity_data = [mean(drd1_data,1); mean(drd2_data,1); ...
                 mean(drd1_mean_per_animal,1); mean(drd2_mean_per_animal,1)];
save(fullfile(drd1_output_dir, 'Activity.mat'), 'activity_data');
save(fullfile(drd1_output_dir, 'Activity_individual.mat'), 'drd1_mean_per_animal');
save(fullfile(drd2_output_dir, 'Activity_individual.mat'), 'drd2_mean_per_animal');

%% 4. Preprocessing
drd1_data(all(isnan(drd1_data) | drd1_data == 0, 2), :) = [];
drd2_data(all(isnan(drd2_data) | drd2_data == 0, 2), :) = [];
save(fullfile(drd1_output_dir, 'DRD1.mat'), 'drd1_data');
save(fullfile(drd2_output_dir, 'DRD2.mat'), 'drd2_data');

%% 5. Normalize and Sort
drd1_z = normalize_data(drd1_data);
drd2_z = normalize_data(drd2_data);
[drd1_sorted, ~] = sort_by_change(drd1_z);
[drd2_sorted, ~] = sort_by_change(drd2_z);

%% 6. Clustering
[drd1_clusters, drd1_means, drd1_labels] = cluster_data(drd1_z, num_clusters, drd1_output_dir, 'DRD1');
[drd2_clusters, drd2_means, drd2_labels] = cluster_data(drd2_z, num_clusters, drd2_output_dir, 'DRD2');

%% 7. Visualization
visualize_data(drd1_sorted, drd1_means, drd1_output_dir, 'DRD1');
visualize_data(drd2_sorted, drd2_means, drd2_output_dir, 'DRD2');

disp('✅ Processing complete.');

%% --- Functions ---

function z = normalize_data(d)
    mu = mean(d(:,1:120), 2);
    sigma = std(d(:,1:120), 0, 2);
    z = (d - mu) ./ sigma;
    z(isnan(z) | isinf(z)) = 0;
end

function [sorted_data, idx] = sort_by_change(z)
    early = mean(z(:, 1:120), 2);
    late = mean(z(:, 160:420), 2);
    [~, idx] = sort(late - early, 'descend');
    sorted_data = z(idx, :);
end

function [clusters, means, labels] = cluster_data(z, k, out_dir, tag)
    change = mean(z(:,160:420), 2) - mean(z(:,1:120),2);
    frequency = sum(z > 2, 2);
    features = [change, frequency];

    [~, score] = pca(features);
    pca_data = score(:,1:2); % 2D for stability

    [labels, ~] = kmeans(pca_data, k, 'Replicates', 10, 'Distance', 'cityblock');
    save(fullfile(out_dir, [tag, '_cluster_labels.mat']), 'labels');

    % Visuals
    figure; gscatter(change, frequency, labels); title([tag ' K-means']);
    saveas(gcf, fullfile(out_dir, [tag '_kmeans.png']));
    close;

    figure; dendrogram(linkage(pca_data, 'ward')); title([tag ' Dendrogram']);
    saveas(gcf, fullfile(out_dir, [tag '_dendrogram.png']));
    close;

    % Elbow
    sse = zeros(1,10);
    for i = 1:10
        [~,~,d] = kmeans(pca_data, i, 'Replicates', 5);
        sse(i) = sum(d);
    end
    figure; plot(1:10, sse, '-o'); title([tag ' Elbow']);
    saveas(gcf, fullfile(out_dir, [tag '_elbow.png']));
    close;

    % Collect clusters
    clusters = cell(k,1);
    for i = 1:k
        clusters{i} = z(labels == i, :);
    end
    means = cell2mat(cellfun(@(x) mean(x,1), clusters, 'UniformOutput', false));
end

function visualize_data(z, cluster_means, out_dir, tag)
    % Cluster mean trace
    figure;
    plot(cluster_means', 'LineWidth', 2);
    xlabel('Time (frames)'); ylabel('Mean Activity');
    title([tag ' Cluster Means']); legend(compose('C%d', 1:size(cluster_means,1)));
    saveas(gcf, fullfile(out_dir, [tag '_cluster_means.png']));
    close;

    % Raw heatmap
    figure;
    imagesc(z); colormap(get_cmap()); caxis([-2 6]); colorbar;
    title([tag ' Heatmap']); xlabel('Time'); ylabel('Neurons');
    saveas(gcf, fullfile(out_dir, [tag '_heatmap.png']));
    close;

    % Blurred heatmap
    blurred = imgaussfilt(z, 2);
    figure;
    imagesc(blurred); colormap(get_cmap()); caxis([-7 7]); colorbar;
    title([tag ' Blurred Heatmap']); xlabel('Time'); ylabel('Neurons');
    saveas(gcf, fullfile(out_dir, [tag '_heatmap_blurred.png']));
    close;
end

function cmap = get_cmap()
    try
        cmap = slanCM('seismic');
    catch
        cmap = colormap('parula'); % fallback
    end
end
