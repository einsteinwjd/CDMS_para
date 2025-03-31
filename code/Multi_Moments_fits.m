clear; close all;
% Moments_fits.m
% This script calculates statistical moments of particle size distributions over time
% from a timetable containing particle number counts in different size bins
path_define;
load([F1_folder,'modeldata_to_timetable.mat']);

% extract the data of first day for test
target_day = 4;
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

% Calculate moments for each time point
for t = 1:numTimePoints
    % Extract particle number distribution at this time
    PN = table2array(simulatedPN(t, 1:end)); % Exclude time column
    
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
figure('Position', [100, 100, 1000, 800]);

for i = 1:length(time_indices)
    t_idx = time_indices(i);
    
    % Original data
    original_data = table2array(simulatedPN(t_idx, 1:end));
    
    % 检测分布是否为多峰分布
    % 计算峰值点
    [peaks, peak_locs] = findpeaks(original_data, 'MinPeakHeight', max(original_data)*0.2, 'MinPeakDistance', 3);
    
    % 使用kurtosis和skewness作为辅助判断
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
    
    % 判断是否需要多峰拟合
    is_multimodal = (length(peaks) > 1) || (kurtosis_value > 3.5);
    
    if is_multimodal
        % 设置初始拟合参数 (假设最多两个模态)
        if length(peaks) >= 2
            % 基于峰位置估计初始参数
            mode1_amp = peaks(1);
            mode1_center = diameterBins(peak_locs(1));
            mode1_width = 0.4; % 初始宽度估计
            
            mode2_amp = peaks(2);
            mode2_center = diameterBins(peak_locs(2));
            mode2_width = 0.4;
            total_number = moments(t_idx, 1);
            initial_params = [mode1_amp, mode1_center, mode1_width, mode2_amp, mode2_center, mode2_width];
        else
            % 如果峰值检测不理想但分布仍然表现为多峰，使用启发式方法设置初始值
            total_number = moments(t_idx, 1);
            initial_params = [total_number*0.6, D_mean(t_idx)*0.7, 0.4, total_number*0.4, D_mean(t_idx)*1.5, 0.5];
        end
        
        % 定义双模态对数正态分布函数
        bimodal_lognormal = @(params, D) (params(1) / (sqrt(2*pi) * params(3) * D)) .* ...
            exp(-(log(D) - log(params(2))).^2 ./ (2 * params(3)^2)) + ...
            (params(4) / (sqrt(2*pi) * params(6) * D)) .* ...
            exp(-(log(D) - log(params(5))).^2 ./ (2 * params(6)^2));
        
        % 定义误差函数
        error_func = @(params) sum((bimodal_lognormal(params, diameterBins) - original_data).^2);
        
        % 非线性优化求解
        options = optimset('Display', 'off', 'MaxIter', 1000);
        
        % 设置参数下界和上界
        lb = [0, minDiameter, 0.1, 0, minDiameter, 0.1]; % 下界
        ub = [total_number*2, maxDiameter, 1.5, total_number*2, maxDiameter, 1.5]; % 上界
        
        [fitted_params, ~] = fmincon(error_func, initial_params, [], [], [], [], lb, ub, [], options);
        
        % 使用拟合参数生成双模态分布
        n_bimodal = bimodal_lognormal(fitted_params, plot_diameters);
        
        % 绘制单峰拟合
        semilogx(plot_diameters, n_lognormal, 'b-', 'LineWidth', 1.5);
        
        % 绘制总的双峰拟合
        semilogx(plot_diameters, n_bimodal, 'r-', 'LineWidth', 2);
        
        % 绘制各个模态的贡献
        mode1 = (fitted_params(1) / (sqrt(2*pi) * fitted_params(3) * plot_diameters)) .* ...
            exp(-(log(plot_diameters) - log(fitted_params(2))).^2 ./ (2 * fitted_params(3)^2));
        mode2 = (fitted_params(4) / (sqrt(2*pi) * fitted_params(6) * plot_diameters)) .* ...
            exp(-(log(plot_diameters) - log(fitted_params(5))).^2 ./ (2 * fitted_params(6)^2));
        
        semilogx(plot_diameters, mode1, 'm--', 'LineWidth', 1);
        semilogx(plot_diameters, mode2, 'g--', 'LineWidth', 1);
        
        legend('原始数据', '单峰拟合', '双峰拟合', '模态1', '模态2');
        
        % 显示拟合参数
        param_text = sprintf('Mode 1: N=%.2e, D_g=%.1f nm, σ_g=%.2f\nMode 2: N=%.2e, D_g=%.1f nm, σ_g=%.2f', ...
            fitted_params(1), fitted_params(2), fitted_params(3), ...
            fitted_params(4), fitted_params(5), fitted_params(6));
        
        % 计算各模态对总浓度的贡献比例
        total_bimodal = fitted_params(1) + fitted_params(4);
        mode1_fraction = fitted_params(1) / total_bimodal * 100;
        mode2_fraction = fitted_params(4) / total_bimodal * 100;
        
        fraction_text = sprintf('Mode 1: %.1f%%, Mode 2: %.1f%%', mode1_fraction, mode2_fraction);
    else
        % 仅单峰分布
        semilogx(plot_diameters, n_lognormal, 'b-', 'LineWidth', 2);
        legend('原始数据', '单峰拟合');
        
        param_text = sprintf('D_g=%.1f nm, σ_g=%.2f', D_g, sigma_g);
        fraction_text = '';
    end
    
    title(['时间点 = ' datestr(timeVector(t_idx), 'HH:MM:SS') ...
           ' (峰值个数=' num2str(length(peaks)) ...
           ', 峭度=' num2str(kurtosis_value, '%.2f') ...
           ', 偏度=' num2str(skewness_value, '%.2f') ')']);
    
    % 添加拟合参数文本
    text(0.05, 0.95, param_text, 'Units', 'normalized', 'FontSize', 8, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'BackgroundColor', [1 1 1 0.7]);
    
    if ~isempty(fraction_text)
        text(0.05, 0.85, fraction_text, 'Units', 'normalized', 'FontSize', 8, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'BackgroundColor', [1 1 1 0.7]);
    end
    
    xlabel('粒径 (nm)');
    ylabel('dN/dlogD');
    grid on;
end

sgtitle('高阶矩分析与多峰粒径分布拟合');

%% 重构随时间变化的热图，考虑多模态分布
figure('Position', [100, 100, 1200, 600]);

% 原始数据热图
subplot(1,2,1);
particleCounts = table2array(simulatedPN(:, 1:numBins));
[Y, X] = meshgrid(diameterBins, 1:numTimePoints);
surf(X, Y, particleCounts);
set(gca, 'YScale', 'log');  % Y轴改为log刻度
view(2); % 平面视图
shading interp;
colorbar;
ylabel('粒径 (nm)');
xlabel('时间点索引');
title('原始粒径分布随时间变化');

% 重构数据热图（考虑多峰分布）
subplot(1,2,2);
reconstructedDistributions = zeros(numTimePoints, length(plot_diameters));

% 对每个时间点检测并拟合分布
for t = 1:numTimePoints
    % 原始数据
    original_data = table2array(simulatedPN(t, 1:end));
    
    % 检测是否为多峰分布
    [peaks, peak_locs] = findpeaks(original_data, 'MinPeakHeight', max(original_data)*0.2, 'MinPeakDistance', 3);
    kurtosis_value = moments(t, 5) * moments(t, 1) / (moments(t, 3)^2);
    
    is_multimodal = (length(peaks) > 1) || (kurtosis_value > 3.5);
    
    if is_multimodal && length(peaks) >= 2
        % 多峰拟合
        % 使用前面类似的方法设置初始参数
        mode1_amp = peaks(1);
        mode1_center = diameterBins(peak_locs(1));
        mode1_width = 0.4;
        
        mode2_amp = peaks(2);
        mode2_center = diameterBins(peak_locs(2));
        mode2_width = 0.4;
        
        initial_params = [mode1_amp, mode1_center, mode1_width, mode2_amp, mode2_center, mode2_width];
        
        % 定义双模态对数正态分布函数
        bimodal_lognormal = @(params, D) (params(1) / (sqrt(2*pi) * params(3) * D)) .* ...
            exp(-(log(D) - log(params(2))).^2 ./ (2 * params(3)^2)) + ...
            (params(4) / (sqrt(2*pi) * params(6) * D)) .* ...
            exp(-(log(D) - log(params(5))).^2 ./ (2 * params(6)^2));
        
        % 定义误差函数
        error_func = @(params) sum((bimodal_lognormal(params, diameterBins) - original_data).^2);
        
        % 优化参数
        options = optimset('Display', 'off', 'MaxIter', 1000);
        total_number = moments(t, 1);
        lb = [0, minDiameter, 0.1, 0, minDiameter, 0.1];
        ub = [total_number*2, maxDiameter, 1.5, total_number*2, maxDiameter, 1.5];
        
        try
            [fitted_params, ~] = fmincon(error_func, initial_params, [], [], [], [], lb, ub, [], options);
            % 使用拟合参数生成重构分布
            reconstructedDistributions(t, :) = bimodal_lognormal(fitted_params, plot_diameters);
        catch
            % 如果优化失败，回退到单峰
            D_g = exp(log(D_mean(t)) - 0.5 * log(geo_std_dev(t)^2));
            sigma_g = geo_std_dev(t);
            
            for j = 1:length(plot_diameters)
                D = plot_diameters(j);
                reconstructedDistributions(t, j) = (moments(t, 1) / (sqrt(2*pi) * log(sigma_g) * D)) * ...
                    exp(-(log(D) - log(D_g))^2 / (2 * log(sigma_g)^2));
            end
        end
    else
        % 单峰拟合
        D_g = exp(log(D_mean(t)) - 0.5 * log(geo_std_dev(t)^2));
        sigma_g = geo_std_dev(t);
        
        for j = 1:length(plot_diameters)
            D = plot_diameters(j);
            reconstructedDistributions(t, j) = (moments(t, 1) / (sqrt(2*pi) * log(sigma_g) * D)) * ...
                exp(-(log(D) - log(D_g))^2 / (2 * log(sigma_g)^2));
        end
    end
end

[Y_recon, X_recon] = meshgrid(plot_diameters, 1:numTimePoints);
surf(X_recon, Y_recon, reconstructedDistributions);
set(gca, 'YScale', 'log');
view(2);
shading interp;
colorbar;
ylabel('粒径 (nm)');
xlabel('时间点索引');
title('基于多峰拟合重构的粒径分布随时间变化');

% 调整颜色映射
colormap(jet);

% 调整颜色范围
for i = [1, 2]
    subplot(1, 2, i);
    if i == 1
        upperLimit = prctile(particleCounts(:), 95);
        caxis([0, upperLimit]);
    else
        upperLimit = prctile(reconstructedDistributions(:), 95);
        caxis([0, upperLimit]); 
    end
    
    c = colorbar;
    c.Label.String = 'dN/dlogD';
end

% 如果时间点适当，添加时间标签
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

% 保存多峰拟合结果
save([F2_folder,'particle_distribution_multimodal_fits.mat'], 'reconstructedDistributions', ...
     'plot_diameters', 'timeVector');

fprintf('多峰分布拟合分析完成并保存到 particle_distribution_multimodal_fits.mat\n');