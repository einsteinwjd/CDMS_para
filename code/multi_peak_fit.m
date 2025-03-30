clear; close all;
% Moments_fits.m
% This script calculates statistical moments of particle size distributions over time
% from a timetable containing particle number counts in different size bins
path_define;
load([F1_folder,'modeldata_to_timetable.mat']);

% extract the data of first day for test
simulatedPN = simulatedPN(1:601,:);


%% 调用多峰分布拟合函数
[fitResults, fittedDistributions] = Multi_peak_distribution_fits(simulatedPN, sim_sizebin);

%% 提取原始数据
timePoints = simulatedPN.Time;
particleCounts = simulatedPN{:, 1:end};
numTimePoints = height(simulatedPN);
numSizeBins = width(particleCounts);

% 创建粒径轴（假设对数等间距）
diameterMin = 1; % nm
diameterMax = 1000; % nm
diameters = logspace(log10(diameterMin), log10(diameterMax), numSizeBins);

%% 可视化比较原始分布与拟合分布
% 选择几个时间点进行比较展示
timePointsToShow = round(linspace(1, numTimePoints, min(6, numTimePoints)));

figure('Position', [100, 100, 1000, 800]);
for i = 1:length(timePointsToShow)
    tp = timePointsToShow(i);
    
    subplot(2, 3, i);
    semilogx(diameters, particleCounts(tp,:), 'bo-', 'DisplayName', '原始数据');
    hold on;
    semilogx(diameters, fittedDistributions(tp,:), 'r-', 'LineWidth', 2, 'DisplayName', '拟合结果');
    
    title(sprintf('时间 = %s', timePoints(tp)));
    xlabel('粒径 (nm)');
    ylabel('数浓度');
    legend('Location', 'best');
    grid on;
end
sgtitle('原始气溶胶粒径分布与多峰拟合比较');

%% 创建随时间变化的热图比较
figure('Position', [100, 100, 1200, 600]);

% 原始数据热图
subplot(2,1,1);
[Y,X] = meshgrid(diameters, 1:numTimePoints);  % 调换X和Y的顺序
surf(X, Y, particleCounts);
set(gca, 'YScale', 'log');  % Y轴改为log刻度，因为现在Y轴是粒径
view(2); % 平面视图
shading interp;
colorbar;
ylabel('粒径 (nm)');  % 修改坐标轴标签
xlabel('时间点索引');
title('原始粒径分布随时间变化');

% 拟合数据热图
subplot(2,1,2);
surf(X, Y, fittedDistributions);
set(gca, 'YScale', 'log');  % Y轴改为log刻度
view(2); % 平面视图
shading interp;
colorbar;
ylabel('粒径 (nm)');  % 修改坐标轴标签
xlabel('时间点索引');
title('拟合粒径分布随时间变化');

%% 分析峰值参数随时间的变化
figure('Position', [100, 100, 1200, 800]);

% 获取所有时间点的峰值数量
maxPeaks = 0;
for t = 1:numTimePoints
    maxPeaks = max(maxPeaks, length(fitResults(t).peaks));
end

% 为每个峰值准备颜色
colors = lines(maxPeaks);

% 峰值中心随时间变化
subplot(3,1,1);
hold on;
for p = 1:maxPeaks
    peakCenters = nan(numTimePoints, 1);
    for t = 1:numTimePoints
        if p <= length(fitResults(t).peaks)
            peakCenters(t) = fitResults(t).peaks(p).center;
        end
    end
    plot(timePoints, peakCenters, 'o-', 'Color', colors(p,:), 'DisplayName', sprintf('峰值 %d', p));
end
xlabel('时间');
ylabel('峰值中心 (nm)');
title('峰值中心随时间演变');
set(gca, 'YScale', 'log');
legend('Location', 'best');
grid on;

% 峰值振幅随时间变化
subplot(3,1,2);
hold on;
for p = 1:maxPeaks
    peakAmplitudes = nan(numTimePoints, 1);
    for t = 1:numTimePoints
        if p <= length(fitResults(t).peaks)
            peakAmplitudes(t) = fitResults(t).peaks(p).amplitude;
        end
    end
    plot(timePoints, peakAmplitudes, 'o-', 'Color', colors(p,:), 'DisplayName', sprintf('峰值 %d', p));
end
xlabel('时间');
ylabel('峰值振幅');
title('峰值振幅随时间演变');
legend('Location', 'best');
grid on;

% 峰值宽度随时间变化
subplot(3,1,3);
hold on;
for p = 1:maxPeaks
    peakWidths = nan(numTimePoints, 1);
    for t = 1:numTimePoints
        if p <= length(fitResults(t).peaks)
            peakWidths(t) = fitResults(t).peaks(p).width;
        end
    end
    plot(timePoints, peakWidths, 'o-', 'Color', colors(p,:), 'DisplayName', sprintf('峰值 %d', p));
end
xlabel('时间');
ylabel('峰值宽度参数');
title('峰值宽度随时间演变');
legend('Location', 'best');
grid on;
