clear; close all;
% Moments_fits.m
% This script calculates statistical moments of particle size distributions over time
% from a timetable containing particle number counts in different size bins
path_define;
load([F1_folder,'modeldata_to_timetable.mat']);

% extract the data of first day for test
target_day = 3;
simulatedPN = simulatedPN(601*(target_day-1)+1:601*target_day,:);

% Check if data exists in workspace, otherwise load it
if ~exist('simulatedPN', 'var')
    % Prompt user to select the file containing the data
    [file, path] = uigetfile('*.mat', 'Select the MAT file containing simulatedPN');
    if file == 0
        error('No file selected');
    end
    load(fullfile(path, file), 'simulatedPN');
end

% Extract time vector
timeVector = simulatedPN.Time;
numTimePoints = height(simulatedPN);

% Extract bin centers (geometric mean diameters)
diameterBins = sim_sizebin;
numBins = length(diameterBins);
minDiameter = min(diameterBins);
maxDiameter = max(diameterBins);

% Generate diameter bin centers (in log space)
diameterBins = logspace(log10(minDiameter), log10(maxDiameter), numBins);

% Initialize arrays to store moments
moments = zeros(numTimePoints, 7); % 0th through 6th moments


% Add a filter to ignore particles below 10nm after the halfway point in time
halfway_point = floor(numTimePoints+1); % filter does not work.

for t = 1:numTimePoints
    % Extract particle number distribution at this time
    PN = table2array(simulatedPN(t, 1:end)); % Exclude time column
    
    % After halfway point, filter out particles below 10nm
    if t > halfway_point
        filter_mask = diameterBins >= 10;
        filtered_PN = PN;
        filtered_PN(~filter_mask) = 0; % Set counts for particles < 10nm to zero
        
        % Calculate moments using filtered data
        for k = 0:6
            moments(t, k+1) = sum(filtered_PN .* diameterBins.^k);
        end
    else
        % Calculate moments as before for first half
        for k = 0:6
            moments(t, k+1) = sum(PN .* diameterBins.^k);
        end
    end
end
% Create a timetable for the moments
momentNames = {'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6'};
momentsTimetable = timetable(timeVector, moments(:,1), moments(:,2), moments(:,3), moments(:,4), ...
                          moments(:,5), moments(:,6), moments(:,7), 'VariableNames', momentNames);

% Calculate derived parameters
% Total number concentration (0th moment)
N_total = moments(:,1);

% Mean diameter (1st moment / 0th moment)
D_mean = moments(:,2) ./ moments(:,1);

% Surface area (proportional to 2nd moment)
S_prop = pi * moments(:,3);

% Volume (proportional to 3rd moment)
V_prop = (pi/6) * moments(:,4);

% Effective diameter (3rd moment / 2nd moment)
D_eff = moments(:,4) ./ moments(:,3);

% Calculate geometric standard deviation
% (Using relationship between moments)
geo_std_dev = sqrt(moments(:,3) .* moments(:,1) ./ (moments(:,2).^2));

% Create plots to visualize the results
figure('Position', [100, 100, 1200, 800]);

% Plot total number concentration
subplot(2,3,1);
plot(timeVector, N_total);
title('Total Number Concentration');
xlabel('Time');
ylabel('N_{total} (count)');
grid on;

% Plot mean diameter
subplot(2,3,2);
plot(timeVector, D_mean);
title('Mean Diameter');
xlabel('Time');
ylabel('D_{mean} (nm)');
grid on;

% Plot surface area
subplot(2,3,3);
plot(timeVector, S_prop);
title('Surface Area (proportional)');
xlabel('Time');
ylabel('Surface Area (nm^2)');
grid on;

% Plot volume
subplot(2,3,4);
plot(timeVector, V_prop);
title('Volume (proportional)');
xlabel('Time');
ylabel('Volume (nm^3)');
grid on;

% Plot effective diameter
subplot(2,3,5);
plot(timeVector, D_eff);
title('Effective Diameter');
xlabel('Time');
ylabel('D_{eff} (nm)');
grid on;

% Plot geometric standard deviation
subplot(2,3,6);
plot(timeVector, geo_std_dev);
title('Geometric Standard Deviation');
xlabel('Time');
ylabel('\sigma_g');
grid on;

% Save results
save([F2_folder,'particle_distribution_moments.mat'], 'momentsTimetable', 'N_total', 'D_mean', ...
     'S_prop', 'V_prop', 'D_eff', 'geo_std_dev');

fprintf('Moment analysis complete. Results saved to particle_distribution_moments.mat\n');

%% check the results by plotting aerosol size distribution by using the moments
% Define a range of diameters for plotting
plot_diameters = logspace(log10(minDiameter), log10(maxDiameter), 100);

% Select a few time points to check
time_indices = round(linspace(1, numTimePoints, 5));

figure('Position', [100, 100, 1000, 600]);

% For each selected time point, reconstruct and plot the distribution
for i = 1:length(time_indices)
    t_idx = time_indices(i);
    
    % Extract moments for this time point
    M0 = moments(t_idx, 1);
    M1 = moments(t_idx, 2);
    M2 = moments(t_idx, 3);
    
    % Assuming lognormal distribution
    % Calculate parameters from moments
    D_g = exp(log(D_mean(t_idx)) - 0.5 * log(geo_std_dev(t_idx)^2));
    sigma_g = geo_std_dev(t_idx);
    
    % Reconstruct the distribution using lognormal equation
    n_reconstructed = zeros(size(plot_diameters));
    for j = 1:length(plot_diameters)
        D = plot_diameters(j);
        n_reconstructed(j) = (M0(1) / (sqrt(2*pi) * log(sigma_g) * D)) * ...
            exp(-(log(D) - log(D_g))^2 / (2 * log(sigma_g)^2));
    end
    
    % Plot reconstructed distribution
    subplot(length(time_indices), 1, i);
    
    % Plot the reconstructed distribution
    semilogx(plot_diameters, n_reconstructed, 'b-', 'LineWidth', 2);
    hold on;
    
    % Plot the original data points for comparison
    semilogx(diameterBins, table2array(simulatedPN(t_idx, 1:end)), 'ro', 'MarkerSize', 6);
    
    title(['Time = ' datestr(timeVector(t_idx), 'HH:MM:SS') ' - Reconstructed vs Original Distribution']);
    xlabel('Diameter (nm)');
    ylabel('dN/dlogD');
    legend('Reconstructed from Moments', 'Original Data');
    grid on;
    hold off;
end

sgtitle('Verification of Moment Analysis by Distribution Reconstruction');

%% Fit multimodal distributions
figure('Position', [100, 100, 1000, 600]);

for i = 1:length(time_indices)
    t_idx = time_indices(i);
    
    % Original data
    original_data = table2array(simulatedPN(t_idx, 1:end));
    
    % Try to fit bimodal lognormal distribution
    % Use higher order moments to characterize distribution shape
    kurtosis_value = moments(t_idx, 5) * moments(t_idx, 1) / (moments(t_idx, 3)^2);
    skewness_value = moments(t_idx, 4) * moments(t_idx, 1) / (moments(t_idx, 3) * moments(t_idx, 2));
    
    subplot(length(time_indices), 1, i);
    semilogx(diameterBins, original_data, 'ko', 'MarkerSize', 6);
    hold on;
    
    % Lognormal fit (single mode)
    D_g = exp(log(D_mean(t_idx)) - 0.5 * log(geo_std_dev(t_idx)^2));
    sigma_g = geo_std_dev(t_idx);
    
    n_lognormal = zeros(size(plot_diameters));
    for j = 1:length(plot_diameters)
        D = plot_diameters(j);
        n_lognormal(j) = (moments(t_idx, 1) / (sqrt(2*pi) * log(sigma_g) * D)) * ...
            exp(-(log(D) - log(D_g))^2 / (2 * log(sigma_g)^2));
    end
    
    semilogx(plot_diameters, n_lognormal, 'b-', 'LineWidth', 1.5);
    
    % If kurtosis and skewness indicate multimodal distribution, add multimodal fit
    if kurtosis_value > 4.2 % Example threshold, should be adjusted based on actual data
        % Implement bimodal lognormal distribution fitting here
        % Need to optimize multiple parameters
        % E.g., mode fractions, geometric mean diameters, and geometric standard deviations
        
        % n_bimodal = mixture_model(plot_diameters, moments(t_idx, :));
        % semilogx(plot_diameters, n_bimodal, 'r-', 'LineWidth', 1.5);
        % legend('Original Data', 'Single Mode Fit', 'Bimodal Fit');
    else
        legend('Original Data', 'Single Mode Fit');
    end
    
    title(['Time = ' datestr(timeVector(t_idx)) ...
           ' (Kurtosis=' num2str(kurtosis_value, '%.2f') ...
           ', Skewness=' num2str(skewness_value, '%.2f') ')']);
    xlabel('Diameter (nm)');
    ylabel('dN/dlogD');
    grid on;
end

sgtitle('Higher-Order Moments Analysis and Fitting of Particle Size Distributions');

%% 创建随时间变化的热图
% 重新组织粒子数据用于热图绘制
particleCounts = table2array(simulatedPN(:, 1:numBins));

% 创建时间与粒径的网格
[Y, X] = meshgrid(diameterBins, 1:numTimePoints);

figure('Position', [100, 100, 1200, 500]);

% 原始数据热图
subplot(1,2,1);
surf(X, Y, particleCounts);
set(gca, 'YScale', 'log');  % Y轴改为log刻度
view(2); % 平面视图
shading interp;
colorbar;
ylabel('粒径 (nm)');
xlabel('时间点索引');
title('粒径分布随时间变化');

% 重构数据热图
subplot(1,2,2);
reconstructedDistributions = zeros(numTimePoints, length(plot_diameters));
for t = 1:numTimePoints
    % 计算重构的分布
    D_g = exp(log(D_mean(t)) - 0.5 * log(geo_std_dev(t)^2));
    sigma_g = geo_std_dev(t);
    
    for j = 1:length(plot_diameters)
        D = plot_diameters(j);
        reconstructedDistributions(t, j) = (moments(t, 1) / (sqrt(2*pi) * log(sigma_g) * D)) * ...
            exp(-(log(D) - log(D_g))^2 / (2 * log(sigma_g)^2));
    end
end

[Y_recon, X_recon] = meshgrid(plot_diameters, 1:numTimePoints);
surf(X_recon, Y_recon, reconstructedDistributions);
set(gca, 'YScale', 'log');  % Y轴改为log刻度
view(2); % 平面视图
shading interp;
colorbar;
ylabel('粒径 (nm)');
xlabel('时间点索引');
title('基于矩重构的粒径分布随时间变化');

% 添加时间标签
colormap(jet); % 使用 jet 色谱以增强对比度

% 调整颜色范围，使热图更清晰
for i = [1, 2]
    subplot(1, 2, i);
    if i == 1
        % 原始数据使用第95百分位数作为上限
        upperLimit = prctile(particleCounts(:), 95);
        caxis([0, upperLimit]);
    else
        % 重构数据使用第95百分位数作为上限
        upperLimit = prctile(reconstructedDistributions(:), 95);
        caxis([0, upperLimit]); 
    end
    
    % 添加颜色条标题
    c = colorbar;
    c.Label.String = 'dN/dlogD';
end

% 添加控制选项，可以选择性地使用对数颜色映射
% 取消下面注释以启用对数颜色尺度
% for i = [1, 2]
%     subplot(1, 2, i);
%     set(gca, 'ColorScale', 'log');
%     caxis([1, max(caxis)]); % 对数尺度下避免0值
% end

% 如果可能，使用实际时间点作为X轴标签
if numTimePoints <= 10
    subplot(1, 2, 1);
    xticks(1:numTimePoints);
    xticklabels(cellstr(datestr(timeVector, 'HH:MM')));
    xtickangle(45);
    
    subplot(1, 2, 2);
    xticks(1:numTimePoints);
    xticklabels(cellstr(datestr(timeVector, 'HH:MM')));
    xtickangle(45);
end


%% J30 comparison between original method and moments approximation
fprintf('Calculating J30 using original and moments-approximated distributions...\n');

% Calculate J30 from original data

    J30_original = Jx_cal(simulatedPN,diameterBins,30);


PN_reconstructed = timetable(simulatedPN.Time,reconstructedDistributions);
    % Calculate J30 from reconstructed distribution
    J30_moments = Jx_cal(PN_reconstructed,diameterBins,30);


% Plot comparison
figure('Position', [100, 100, 1000, 500]);
plot(timeVector, J30_original, 'b-', 'LineWidth', 2);
hold on;
plot(timeVector, J30_moments, 'r--', 'LineWidth', 2);
hold off;
title('J30 Comparison: Original vs. Moments Approximation');
xlabel('Time');
ylabel('J30 (cm^{-3} s^{-1})');
legend('Original Method', 'Moments Approximation');
grid on;

% Calculate error metrics
relative_error = abs(J30_moments - J30_original) ./ J30_original * 100;
mean_rel_error = mean(relative_error, 'omitnan');
max_rel_error = max(relative_error, [], 'omitnan');

% Add text annotation with error metrics
text(0.05, 0.9, sprintf('Mean Relative Error: %.2f%%\nMax Relative Error: %.2f%%', ...
     mean_rel_error, max_rel_error), 'Units', 'normalized');

% Save results for further analysis
save([F2_folder,'J30_comparison.mat'], 'J30_original', 'J30_moments', ...
     'mean_rel_error', 'max_rel_error');

fprintf('J30 comparison complete. Results saved to J30_comparison.mat\n');