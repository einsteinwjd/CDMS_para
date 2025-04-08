% clear; close all;
function [] = Moments_fits_obs(target_day,path)
close all;
% Moments_fits.m
% This script calculates statistical moments of particle size distributions over time
% from a timetable containing particle number counts in different size bins
eval(path);
load([F1_folder,'timetable_data_nj.mat']);

% extract the data of first day for test
% target_day = 1;
obsPN = tt_psd(120*(target_day-1)+1:120*target_day,:);

% Check if data exists in workspace, otherwise load it
if ~exist('obsPN', 'var')
    % Prompt user to select the file containing the data
    [file, path] = uigetfile('*.mat', 'Select the MAT file containing simulatedPN');
    if file == 0
        error('No file selected');
    end
    load(fullfile(path, file), 'obsPN');
end

% Extract time vector
timeVector = obsPN.Time;
numTimePoints = height(obsPN);

% Extract bin centers from column names of simulatedPN
varNames = obsPN.Properties.VariableNames;
diameterBins = zeros(1, length(varNames));

for i = 1:length(varNames)
    % Extract numeric part after 'd' character (e.g., 'd1.4' -> 1.4)
    varName = varNames{i};
    if startsWith(varName, 'd')
        diameterStr = extractAfter(varName, 'd');
        diameterBins(i) = str2double(diameterStr);
    end
end

% Remove any zero values (in case there were non-diameter columns)
diameterBins = diameterBins(diameterBins > 0);
numBins = length(diameterBins);
minDiameter = min(diameterBins);
maxDiameter = max(diameterBins);

% Generate diameter bin centers (in log space)
diameterBins = logspace(log10(minDiameter), log10(maxDiameter), numBins);

% Initialize arrays to store moments
moments = zeros(numTimePoints, 7); % 0th through 6th moments

% Calculate moments for each time point
for t = 1:numTimePoints
    % Extract particle number distribution at this time
    PN = table2array(obsPN(t, 1:end)); % Exclude time column
    
    % Calculate moments (k = 0, 1, 2, 3, 4, 5, 6)
    for k = 0:6
        moments(t, k+1) = sum(PN .* diameterBins.^k);
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
save([F2_folder,'particle_distribution_moments_obs.mat'], 'momentsTimetable', 'N_total', 'D_mean', ...
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
    semilogx(diameterBins, table2array(obsPN(t_idx, 1:end)), 'ro', 'MarkerSize', 6);
    
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
    original_data = table2array(obsPN(t_idx, 1:end));
    
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
particleCounts = table2array(obsPN(:, 1:numBins));

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
    tick_index = 1:24:numTimePoints;
    subplot(1, 2, 1);
    % xlim([1,numTimePoints+1]);
    xticks(tick_index);
    xticklabels(cellstr(datestr(timeVector(tick_index), 'HH:MM')));
    xtickangle(45);

    subplot(1, 2, 2);
    % xlim([1,1201]);
    xticks(tick_index);
    xticklabels(cellstr(datestr(timeVector(tick_index), 'HH:MM')));
    xtickangle(45);

%% J30 comparison between original method and moments approximation
fprintf('计算J30的观测值和moments近似值时间序列对比...\n');

% 计算原始数据的J30
J30_original = Jx_cal(obsPN, diameterBins, 30);

% 创建基于moments重构的粒径分布的timetable
PN_reconstructed = array2timetable(reconstructedDistributions, 'RowTimes', timeVector);

% 计算重构分布的J30
J30_moments = Jx_cal(PN_reconstructed, plot_diameters, 30);

% 绘制对比图
figure('Position', [100, 100, 1000, 500]);
plot(timeVector, J30_original, 'b-', 'LineWidth', 2);
hold on;
plot(timeVector, J30_moments, 'r--', 'LineWidth', 2);
hold off;
title('J30对比: 观测值 vs. Moments近似值');
xlabel('时间');
ylabel('J30 (cm^{-3} s^{-1})');
legend('观测值', 'Moments近似值');
grid on;
datetick('x', 'HH:MM', 'keepticks');

% 计算误差指标
valid_indices = (J30_original > 0) & ~isnan(J30_original) & ~isnan(J30_moments);
relative_error = abs(J30_moments(valid_indices) - J30_original(valid_indices)) ./ J30_original(valid_indices) * 100;
mean_rel_error = mean(relative_error, 'omitnan');
max_rel_error = max(relative_error, [], 'omitnan');

% 添加误差指标文本注释
text(0.05, 0.9, sprintf('平均相对误差: %.2f%%\n最大相对误差: %.2f%%', ...
     mean_rel_error, max_rel_error), 'Units', 'normalized');

% 保存图片到Fig_folder (300 DPI)，加入target_day信息
saveas(gcf, [Fig_folder, sprintf('J30_comparison_timeseries_obs_day%d.fig', target_day)]);
print(gcf, [Fig_folder, sprintf('J30_comparison_timeseries_obs_day%d_300dpi', target_day)], '-dpng', '-r300');

% 保存结果供进一步分析
save([F2_folder,'J30_comparison_obs.mat'], 'J30_original', 'J30_moments', ...
     'mean_rel_error', 'max_rel_error', 'timeVector');

fprintf('J30对比完成。结果已保存至J30_comparison.mat\n');

% 额外创建一个绘图展示粒径分布和30nm的位置
figure('Position', [100, 100, 1000, 400]);
subplot(1, 2, 1);
% 选择一个时间点展示
sample_time_idx = round(numTimePoints/2);
semilogx(diameterBins, particleCounts(sample_time_idx, :), 'bo-', 'LineWidth', 1.5);
hold on;
semilogx(plot_diameters, reconstructedDistributions(sample_time_idx, :), 'r-', 'LineWidth', 2);
% 标记30nm位置
plot([30 30], ylim, 'k--', 'LineWidth', 1.5);
xlabel('粒径 (nm)');
ylabel('dN/dlogD');
title(['时间点: ', datestr(timeVector(sample_time_idx), 'HH:MM'), ' 的粒径分布']);
legend('观测分布', 'Moments重构', '30 nm', 'Location', 'best');
grid on;

% 绘制不同时间点的J30通量变化
subplot(1, 2, 2);
times_to_show = round(linspace(1, numTimePoints, min(6, numTimePoints)));
hold on;
for i = 1:length(times_to_show)
    t_idx = times_to_show(i);
    plot(diameterBins, particleCounts(t_idx, :) ./ max(particleCounts(t_idx, :)), 'LineWidth', 1.5);
end
% 标记30nm位置
plot([30 30], ylim, 'k--', 'LineWidth', 1.5);
set(gca, 'XScale', 'log');
title('不同时间点的归一化粒径分布');
xlabel('粒径 (nm)');
ylabel('归一化dN/dlogD');
legend_str = cellstr(datestr(timeVector(times_to_show), 'HH:MM'));
legend([legend_str; '30 nm'], 'Location', 'best');
grid on;

% 保存图片到Fig_folder (300 DPI)
saveas(gcf, [Fig_folder, sprintf('J30_distribution_comparison_obs_day%d.fig', target_day)]);
print(gcf, [Fig_folder, sprintf('J30_distribution_comparison_obs_day%d_300dpi', target_day)], '-dpng', '-r300');

% 同样为之前生成的图保存300 DPI版本
% 保存矩分析主图
figure(1); % 假设这是矩分析的主图
saveas(gcf, [Fig_folder, sprintf('moments_analysis_obs_day%d.fig', target_day)]);
print(gcf, [Fig_folder, sprintf('moments_analysis_obs_day%d_300dpi', target_day)], '-dpng', '-r300');

% 保存分布重构验证图
figure(2); % 假设这是分布重构验证图
saveas(gcf, [Fig_folder, sprintf('distribution_reconstruction_obs_day%d.fig', target_day)]);
print(gcf, [Fig_folder, sprintf('distribution_reconstruction_obs_day%d_300dpi', target_day)], '-dpng', '-r300');

% 保存高阶矩分析图
figure(3); % 假设这是高阶矩分析图
saveas(gcf, [Fig_folder, sprintf('higher_order_moments_analysis_obs_day%d.fig', target_day)]);
print(gcf, [Fig_folder, sprintf('higher_order_moments_analysis_obs_day%d_300dpi', target_day)], '-dpng', '-r300');

% 保存热图
figure(4); % 假设这是热图
saveas(gcf, [Fig_folder, sprintf('particle_distribution_heatmap_obs_day%d.fig', target_day)]);
print(gcf, [Fig_folder, sprintf('particle_distribution_heatmap_obs_day%d_300dpi', target_day)], '-dpng', '-r300');

fprintf('所有图片已保存为300 DPI的PNG格式到 %s，文件名中包含第%d天信息\n', Fig_folder, target_day);

end