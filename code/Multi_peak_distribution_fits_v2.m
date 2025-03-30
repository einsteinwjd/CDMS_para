function [fitResults, fittedDistributions] = Multi_peak_distribution_fits_v2(simulatedPN, sim_sizebin)
% Multi-peak distributino Fits multi-peak distributions to particle size data over time
%
% Inputs:
%   simulatedPN - timetable with time in first column and particle counts in 80 columns
%   sim_sizebin - vector containing geometric mean diameters of size bins
%
% Outputs:
%   fitResults - structure containing fitting parameters for each timepoint
%   fittedDistributions - fitted distributions for each timepoint

% Extract data from timetable
timePoints = simulatedPN.Time;
particleCounts = simulatedPN{:, 1:end};
numTimePoints = height(simulatedPN);
numSizeBins = width(particleCounts);

% Use sim_sizebin for diameter axis instead of creating a new one
diameters = sim_sizebin;

% Verify dimensions match
if length(diameters) ~= numSizeBins
    error('Length of sim_sizebin (%d) must match the number of size bins in simulatedPN (%d)', length(diameters), numSizeBins);
end

% Storage for results
fitResults = struct('time', cell(numTimePoints, 1), 'peaks', cell(numTimePoints, 1), 'gof', cell(numTimePoints, 1));
fittedDistributions = zeros(numTimePoints, numSizeBins);

% Define parallel pool
if isempty(gcp('nocreate'))
    parpool('local'); % Create local parallel pool
end

% Process each time point in parallel
parfor t = 1:numTimePoints
    % Current particle size distribution
    currentDist = particleCounts(t,:);
    
    % Determine number of peaks
    numPeaks = estimateNumberOfPeaks(currentDist, diameters);
    
    % Create fitting function with appropriate number of peaks
    fitType = createMultiPeakFitType(numPeaks);
    
    % Set initial parameters based on distribution
    initialParams = estimateInitialParameters(currentDist, diameters, numPeaks);
    
    % Set parameter bounds - prevent parameters from deviating during fitting
    lower = [];
    upper = [];
    for i = 1:numPeaks
        % a: amplitude lower bound is 0, upper bound is 3 times the maximum value
        % b: log center within the log diameter range
        % c: width limited to a reasonable range
        lower = [lower, 0, log(min(diameters))*0.9, 0.05];
        upper = [upper, max(currentDist)*3, log(max(diameters))*1.1, 2];
    end
    
    % Create local fitoptions object (avoid issues in parfor)
    options = fitoptions('Method', 'NonlinearLeastSquares',...
                    'Lower', lower,...
                    'Upper', upper,...
                    'StartPoint', initialParams,...
                    'MaxIter', 1000,...
                    'TolFun', 1e-8,...
                    'TolX', 1e-8);
    
    % Use simulated annealing-like multiple fitting strategy
    bestGof = Inf;
    bestFit = [];
    
    % Create alternative initial parameter sets
    paramSets = {};
    weightSets = [];
    
    % Add estimated initial parameters
    paramSets{1} = initialParams;
    weightSets(1) = 1.0;
    
    % Create additional parameter sets for overlapping peak regions
    if numPeaks >= 2
        % Detect potential overlapping regions
        potentialOverlaps = findPotentialOverlaps(currentDist, diameters);
        
        if ~isempty(potentialOverlaps)
            % Generate additional parameter sets focusing on overlapping regions
            for i = 1:min(2, length(potentialOverlaps))
                overlapParams = adjustParamsForOverlap(initialParams, numPeaks, potentialOverlaps(i), diameters);
                paramSets{end+1} = overlapParams;
                weightSets(end+1) = 0.8;
            end
        end
    end
    
    % Multiple attempts with different starting points
    maxAttempts = 5 + numPeaks;  % More peaks require more attempts
    
    for attempt = 1:maxAttempts
        try
            % Dynamically select initial parameters
            if attempt <= length(paramSets)
                % Use preset parameter sets
                options.StartPoint = paramSets{attempt};
            else
                % Perturb based on previous best result using annealing strategy
                if ~isempty(bestFit)
                    baseParams = coeffvalues(bestFit);
                    
                    % Increase perturbation range initially, then decrease
                    if attempt < maxAttempts*0.7
                        perturbRange = 0.2 + 0.3*(attempt/maxAttempts);
                    else
                        perturbRange = 0.3 - 0.2*((attempt-maxAttempts*0.7)/(maxAttempts*0.3));
                    end
                    
                    perturbation = 1 + perturbRange * (rand(size(baseParams))*2-1);
                    options.StartPoint = baseParams .* perturbation;
                    
                    % Ensure parameters are within reasonable bounds
                    options.StartPoint = max(min(options.StartPoint, upper*0.95), lower*1.05);
                else
                    % Slightly perturb initial parameters
                    perturbedParams = initialParams .* (0.7 + 0.6*rand(size(initialParams)));
                    options.StartPoint = perturbedParams;
                end
            end
            
            % Adjust optimization settings - gradually increase precision
            if attempt < maxAttempts/2
                options.MaxIter = 800;
                options.TolFun = 1e-7;
                options.TolX = 1e-7;
            else
                options.MaxIter = 1200;
                options.TolFun = 1e-8;
                options.TolX = 1e-8;
            end
            
            % Perform fitting
            [fitResult, gof] = fit(diameters', currentDist', fitType, options);
            
            % Evaluate fit quality - consider RMSE and physical significance of peaks
            qualityScore = evaluateFitQuality(fitResult, gof, currentDist, diameters, numPeaks);
            
            % Save best fit result
            if qualityScore < bestGof
                bestGof = qualityScore;
                bestFit = fitResult;
            end
        catch ME
            % Log error but continue attempts
            warning('Time point %d, attempt %d failed: %s', t, attempt, ME.message);
            continue;
        end
    end
    
    % Use best fit result or attempt to reduce complexity
    if ~isempty(bestFit)
        fitResult = bestFit;
        gof.rmse = bestGof;
    else
        % If all attempts fail, try different strategies
        success = false;
        
        % Strategy 1: Reduce number of peaks
        if numPeaks > 1
            try
                numPeaks = numPeaks - 1;
                fitType = createMultiPeakFitType(numPeaks);
                initialParams = estimateInitialParameters(currentDist, diameters, numPeaks);
                options.StartPoint = initialParams;
                options.Lower = lower(1:3*numPeaks);
                options.Upper = upper(1:3*numPeaks);
                [fitResult, gof] = fit(diameters', currentDist', fitType, options);
                success = true;
            catch
                % Continue to next strategy
            end
        end
        
        % Strategy 2: Use simpler overlapping peak-specific model
        if ~success && numPeaks >= 2
            try
                specialType = createOverlappingPeakFitType(2); % Special overlapping peak model
                specialParams = estimateOverlappingParams(currentDist, diameters);
                options.StartPoint = specialParams;
                [fitResult, gof] = fit(diameters', currentDist', specialType, options);
                success = true;
            catch
                % Continue to next strategy
            end
        end
        
        % Final fallback: Use single-peak Gaussian model
        if ~success
            warning('Time point %d fitting failed, using single-peak Gaussian model', t);
            fitType = fittype('a1*exp(-((log(x)-b1)^2)/(2*c1^2))', 'independent', 'x');
            options.StartPoint = [max(currentDist), log(diameters(round(end/2))), 0.5];
            [fitResult, gof] = fit(diameters', currentDist', fitType, options);
            numPeaks = 1;
        end
    end
    
    % Post-process - check and correct unreasonable fit results
    fitResult = postProcessFit(fitResult, currentDist, diameters, numPeaks);
    
    % Store results (note that in parfor, each element must be stored separately)
    localFitResult = struct('time', timePoints(t), ...
                           'peaks', extractPeakParameters(fitResult, numPeaks), ...
                           'gof', gof);
    
    % Store results, handle each index separately in parfor
    fitResults(t) = localFitResult;
    fittedDistributions(t,:) = fitResult(diameters);
end

% Analyze time evolution of parameters after parallel processing
analyzeParameterEvolution(fitResults, timePoints);

end

function numPeaks = estimateNumberOfPeaks(distribution, diameters)
% 使用多尺度分析方法估计分布中的峰值数量，特别优化对交叠峰的处理
    % 计算对数坐标下的分布 - 有助于分辨接近的峰值
    logDiameters = log10(diameters);
    logSpacing = mean(diff(logDiameters));
    
    % 使用多种平滑宽度进行分析
    smoothWidths = [3, 5, 9];
    allPeaks = [];
    allLocs = [];
    allProms = [];
    
    % 在不同平滑尺度下寻找峰值
    for w = smoothWidths
        smoothDist = smoothdata(distribution, 'gaussian', w);
        [pks, locs, widths, proms] = findpeaks(smoothDist, 'MinPeakDistance', 2);
        
        if ~isempty(pks)
            % 只保留相对突出的峰值
            maxHeight = max(smoothDist);
            validIdx = proms > (maxHeight * 0.03);  % 降低阈值以检测交叠峰
            
            allPeaks = [allPeaks, pks(validIdx)];
            allLocs = [allLocs, locs(validIdx)];
            allProms = [allProms, proms(validIdx)];
        end
    end
    
    % 如果常规方法找不到足够峰值，尝试使用二阶导数信息
    if length(allPeaks) <= 1
        % 二阶导数可以帮助识别交叠峰
        smooth_dist = smoothdata(distribution, 'gaussian', 5);
        d2 = diff(diff(smooth_dist));
        d2 = [0, d2, 0];  % 补齐长度
        
        % 寻找二阶导数的负峰值（对应原信号的峰值）
        [~, potentialPeaks] = findpeaks(-d2, 'MinPeakProminence', max(abs(d2))*0.05);
        
        for i = 1:length(potentialPeaks)
            idx = potentialPeaks(i);
            if idx > 1 && idx < length(distribution) && ~ismember(idx, allLocs)
                % 确认是真实峰值而不只是噪声
                if (distribution(idx) > distribution(idx-1) || distribution(idx) > distribution(idx+1))
                    allPeaks = [allPeaks; distribution(idx)];
                    allLocs = [allLocs; idx];
                    allProms = [allProms; max(distribution(idx)-min(distribution(max(1,idx-3):min(end,idx+3))), 0.1)];
                end
            end
        end
    end
    
    % 聚类合并接近的峰值
    if length(allLocs) > 1
        clusterDist = 3; % 接近程度阈值
        [sortedLocs, idx] = sort(allLocs);
        sortedPeaks = allPeaks(idx);
        sortedProms = allProms(idx);
        
        i = 1;
        mergedLocs = [];
        mergedPeaks = [];
        mergedProms = [];
        
        while i <= length(sortedLocs)
            cluster = sortedLocs(i);
            maxPeak = sortedPeaks(i);
            maxProm = sortedProms(i);
            j = i + 1;
            
            % 查找临近峰值
            while j <= length(sortedLocs) && abs(sortedLocs(j) - sortedLocs(i)) <= clusterDist
                if sortedProms(j) > maxProm
                    maxPeak = sortedPeaks(j);
                    maxProm = sortedProms(j);
                    cluster = sortedLocs(j);
                end
                j = j + 1;
            end
            
            mergedLocs = [mergedLocs; cluster];
            mergedPeaks = [mergedPeaks; maxPeak];
            mergedProms = [mergedProms; maxProm];
            i = j;
        end
        
        allLocs = mergedLocs;
        allPeaks = mergedPeaks;
        allProms = mergedProms;
    end
    
    % 根据突出度排序并选择前几个峰值
    [sortedProms, idx] = sort(allProms, 'descend');
    allLocs = allLocs(idx);
    
    % 确定峰值数量，限制在1-4个之间
    numPeaks = min(4, max(1, min(3, length(allLocs))));
    
    % 调试用的可视化代码
    % figure; plot(diameters, distribution, 'b-'); hold on;
    % sm_dist = smoothdata(distribution, 'gaussian', 5);
    % plot(diameters, sm_dist, 'k--');
    % plot(diameters(allLocs(1:min(numPeaks,length(allLocs)))), 
    %      distribution(allLocs(1:min(numPeaks,length(allLocs)))), 'ro', 'MarkerFaceColor', 'r');
    % set(gca, 'XScale', 'log'); title(['检测到 ', num2str(numPeaks), ' 个峰值']);
end

function fitType = createMultiPeakFitType(numPeaks)
% Create a fit type string for multi-peak lognormal distribution
    expression = '';
    for i = 1:numPeaks
        if i > 1
            expression = [expression ' + '];
        end
        % a*exp(-((log(x)-b)^2)/(2*c^2)) for each peak
        expression = [expression sprintf('a%d*exp(-((log(x)-b%d)^2)/(2*c%d^2))', i, i, i)];
    end
    fitType = fittype(expression, 'independent', 'x');
end

function params = estimateInitialParameters(distribution, diameters, numPeaks)
% 优化的初始参数估计函数，特别加强对交叠峰的处理
    params = [];
    
    % 使用不同平滑宽度进行多尺度分析
    smoothDist1 = smoothdata(distribution, 'gaussian', 3);  % 细节特征
    smoothDist2 = smoothdata(distribution, 'gaussian', 7);  % 主要特征
    
    % 尝试对明显的交叠峰进行预分离
    separatedDist = attemptPeakSeparation(distribution, diameters);
    
    % 在原始分布和平滑分布上寻找峰值
    [pks1, locs1] = findpeaks(smoothDist1, 'MinPeakProminence', max(smoothDist1)*0.04);
    [pks2, locs2] = findpeaks(smoothDist2, 'MinPeakProminence', max(smoothDist2)*0.06);
    [pks3, locs3] = findpeaks(separatedDist, 'MinPeakProminence', max(separatedDist)*0.05);
    
    % 收集所有可能的峰值
    allPeaks = []; allLocs = [];
    
    if ~isempty(locs1), allPeaks = [allPeaks, pks1]; allLocs = [allLocs, locs1]; end
    if ~isempty(locs2), allPeaks = [allPeaks, pks2]; allLocs = [allLocs, locs2]; end
    if ~isempty(locs3), allPeaks = [allPeaks, pks3]; allLocs = [allLocs, locs3]; end
    
    % 使用二阶导数信息找出可能被掩盖的峰值
    if length(allPeaks) < numPeaks || numPeaks >= 2
        d2 = diff(diff(smoothDist1)); d2 = [0, d2, 0];
        [~, d2Peaks] = findpeaks(-d2, 'MinPeakProminence', 0.05*max(abs(d2)));
        
        for i = 1:length(d2Peaks)
            idx = d2Peaks(i);
            if ~ismember(idx, allLocs) && idx > 1 && idx < length(distribution)
                allLocs = [allLocs, idx];
                allPeaks = [allPeaks, distribution(idx)];
            end
        end
    end
    
    % 如果仍然找不到足够峰值，添加基于分布形状的猜测位置
    if length(allPeaks) < numPeaks
        % 识别分布中的"肩部"(可能的隐藏峰)
        shoulderPoints = findShoulderPoints(distribution, diameters);
        for i = 1:length(shoulderPoints)
            if ~ismember(shoulderPoints(i), allLocs)
                allLocs = [allLocs, shoulderPoints(i)];
                allPeaks = [allPeaks, distribution(shoulderPoints(i))];
            end
        end
        
        % 如果依然不够，添加均匀分布的峰值
        if length(allPeaks) < numPeaks
            additionalPoints = round(linspace(10, length(diameters)-10, numPeaks+2));
            additionalPoints = additionalPoints(2:end-1);
            for i = 1:length(additionalPoints)
                if ~ismember(additionalPoints(i), allLocs)
                    allLocs = [allLocs, additionalPoints(i)];
                    allPeaks = [allPeaks, distribution(additionalPoints(i))*0.8];
                end
            end
        end
    end
    
    % 去除重复点并排序
    [allLocs, uniqueIdx] = unique(allLocs);
    allPeaks = allPeaks(uniqueIdx);
    [allPeaks, sortIdx] = sort(allPeaks, 'descend');
    allLocs = allLocs(sortIdx);
    
    % 取前numPeaks个峰
    pks = allPeaks(1:min(numPeaks, length(allPeaks)));
    locs = allLocs(1:min(numPeaks, length(allLocs)));
    
    % 如果峰位数量不足，补充
    if length(locs) < numPeaks
        missingPeaks = numPeaks - length(locs);
        logRange = log(max(diameters)) - log(min(diameters));
        positions = linspace(log(min(diameters))+0.2*logRange, log(max(diameters))-0.2*logRange, missingPeaks+2);
        positions = positions(2:end-1);
        
        for i = 1:missingPeaks
            pks = [pks, max(distribution)*0.3];
            
            % 找最接近positions(i)的位置
            [~, posIdx] = min(abs(log(diameters) - positions(i)));
            locs = [locs, posIdx];
        end
    end
    
    % 为交叠峰优化参数估计
    for i = 1:numPeaks
        if i <= length(pks)
            a = pks(i);  % 高度
            peakX = diameters(locs(i));
            b = log(peakX);  % 对数中心
            
            % 获取考虑邻近峰影响的宽度估计
            c = estimatePeakWidth(distribution, diameters, locs(i), locs, i);
            
            % 对交叠峰区域进行特殊处理
            isOverlapping = false;
            for j = 1:numPeaks
                if j ~= i && abs(locs(j) - locs(i)) < 10
                    isOverlapping = true;
                    break;
                end
            end
            
            if isOverlapping
                % 对交叠峰区域的宽度进行额外修正
                c = c * (0.8 + 0.2*rand());  % 略微扰动以避免局部最小值
                
                % 高度可能受影响，尝试估计真实高度
                a = estimateOverlappingPeakHeight(distribution, diameters, locs, i);
            end
            
            % 确保宽度在合理范围内
            c = min(max(c, 0.1), 1.2);
        else
            % 默认参数
            a = max(distribution) * 0.2;
            b = log(diameters(round(length(diameters)/2)));
            c = 0.4;
        end
        
        params = [params, a, b, c];
    end
end

function separated = attemptPeakSeparation(distribution, diameters)
% 增强交叠峰分离能力
    separated = distribution;
    smoothed = smoothdata(distribution, 'gaussian', 5);
    
    % 检测可能的交叠区域(曲率变化区域)
    d2 = diff(diff(smoothed));
    d2 = [0, d2, 0];
    
    % 使用更强大的交叠峰检测方法
    % 寻找二阶导数的极值(可能是交叠区)
    [~, extremaIdx] = findpeaks(abs(d2), 'MinPeakProminence', max(abs(d2))*0.08);
    inflections = extremaIdx;
    
    % 通过一阶导数验证交叠区域
    d1 = diff(smoothed);
    d1 = [d1(1), d1];
    
    % 寻找导数变化异常的区域
    validInflections = [];
    for i = 1:length(inflections)
        idx = inflections(i);
        if idx > 5 && idx < length(distribution)-5
            % 检查前后区间的一阶导数变化率
            preSlope = mean(d1(idx-5:idx-1));
            postSlope = mean(d1(idx+1:idx+5));
            
            % 如果导数变化异常，更可能是交叠区
            if preSlope * postSlope < 0 || abs(preSlope/postSlope) > 2
                validInflections = [validInflections, idx];
            end
        end
    end
    
    if ~isempty(validInflections)
        inflections = validInflections;
    end
    
    % 在交叠区域应用自适应增强分离
    for i = 1:length(inflections)
        idx = inflections(i);
        if idx > 3 && idx < length(distribution)-3
            % 计算曲率强度，决定分离程度
            curvatureStrength = abs(d2(idx)) / mean(abs(d2));
            
            % 动态决定分离深度 - 较强的曲率需要更强的分离
            separationDepth = min(max(0.75, 0.85 - 0.1*(curvatureStrength-1)), 0.95);
            
            % 应用递减型分离，避免突变
            separated(idx) = separated(idx) * separationDepth;
            separated(idx-1) = separated(idx-1) * (separationDepth + 0.05);
            separated(idx+1) = separated(idx+1) * (separationDepth + 0.05);
            separated(idx-2) = separated(idx-2) * (separationDepth + 0.1);
            separated(idx+2) = separated(idx+2) * (separationDepth + 0.1);
        end
    end
    
    % 平滑处理分离后的分布，避免引入假峰
    separated = smoothdata(separated, 'gaussian', 2);
end

function shoulderPoints = findShoulderPoints(distribution, diameters)
% 识别分布中的"肩部"(可能是被主峰掩盖的次峰)
    shoulderPoints = [];
    smoothDist = smoothdata(distribution, 'gaussian', 5);
    
    % 计算一阶导数
    d1 = diff(smoothDist);
    d1 = [d1(1), d1];
    
    % 寻找导数变化剧烈的区域
    d2 = diff(d1);
    d2 = [d2(1), d2];
    
    % 找出局部曲率特征明显的位置
    for i = 5:length(d2)-5
        if abs(d2(i)) > 3*std(d2) && distribution(i) > 0.1*max(distribution)
            % 确认是肩部而非噪声
            window = i-4:i+4;
            localDist = distribution(window);
            localD2 = d2(window);
            
            if std(localDist) < 0.2*max(distribution) && std(localD2) > std(d2)
                shoulderPoints = [shoulderPoints, i];
            end
        end
    end
end

function c = estimatePeakWidth(distribution, diameters, peakLoc, allLocs, currentPeakIdx)
% 考虑邻近峰影响的宽度估计
    peakHeight = distribution(peakLoc);
    halfMax = peakHeight / 2;
    
    % 确定左右搜索边界(避免受邻近峰影响)
    leftBound = 1;
    rightBound = length(distribution);
    
    % 邻近峰位置会限制搜索范围
    for j = 1:length(allLocs)
        if j ~= currentPeakIdx
            if allLocs(j) < peakLoc && allLocs(j) > leftBound
                midPoint = round((allLocs(j) + peakLoc)/2);
                leftBound = max(leftBound, midPoint);
            elseif allLocs(j) > peakLoc && allLocs(j) < rightBound
                midPoint = round((allLocs(j) + peakLoc)/2);
                rightBound = min(rightBound, midPoint);
            end
        end
    end
    
    % 在有限范围内搜索半高宽
    leftIdx = peakLoc;
    while leftIdx > leftBound && distribution(leftIdx) > halfMax
        leftIdx = leftIdx - 1;
    end
    
    rightIdx = peakLoc;
    while rightIdx < rightBound && distribution(rightIdx) > halfMax
        rightIdx = rightIdx + 1;
    end
    
    % 估计半高宽
    if leftIdx == leftBound && distribution(leftIdx) > halfMax
        leftWidth = 2 * (peakLoc - leftIdx);
    else
        if leftIdx < peakLoc
            leftPos = leftIdx + (halfMax - distribution(leftIdx)) / ...
                (distribution(leftIdx+1) - distribution(leftIdx));
        else
            leftPos = leftIdx;
        end
        leftWidth = peakLoc - leftPos;
    end
    
    if rightIdx == rightBound && distribution(rightIdx) > halfMax
        rightWidth = 2 * (rightIdx - peakLoc);
    else
        if rightIdx > peakLoc && rightIdx <= length(distribution)
            rightPos = (rightIdx-1) + (halfMax - distribution(rightIdx-1)) / ...
                (distribution(rightIdx) - distribution(rightIdx-1));
        else
            rightPos = rightIdx;
        end
        rightWidth = rightPos - peakLoc;
    end
    
    % 计算对数空间中的宽度
    logWidth = (leftWidth + rightWidth) * (log(diameters(2)) - log(diameters(1)));
    c = logWidth / 2.355;  % FWHM → 标准差
end

function realHeight = estimateOverlappingPeakHeight(distribution, diameters, locs, currentPeakIdx)
% 估计交叠峰的真实高度
    peakLoc = locs(currentPeakIdx);
    obsHeight = distribution(peakLoc);
    
    % 默认高度就是观测高度
    realHeight = obsHeight;
    
    % 检查是否处于交叠区域
    isOverlapped = false;
    nearbyPeaks = [];
    
    for j = 1:length(locs)
        if j ~= currentPeakIdx && abs(locs(j) - peakLoc) < 10
            isOverlapped = true;
            nearbyPeaks = [nearbyPeaks, j];
        end
    end
    
    if isOverlapped && ~isempty(nearbyPeaks)
        % 粗略估计：假设每个峰叠加贡献了一部分高度
        % 根据距离估计贡献比例
        contributionFactor = 1.0; % 默认
        
        for j = nearbyPeaks
            dist = abs(locs(j) - peakLoc);
            
            % 距离越近，贡献越大
            if dist <= 3
                contributionFactor = contributionFactor + 0.5;
            elseif dist <= 7
                contributionFactor = contributionFactor + 0.2;
            else
                contributionFactor = contributionFactor + 0.1;
            end
        end
        
        % 调整估计高度
        realHeight = obsHeight * (1 + 0.2*(contributionFactor-1));
    end
end

function potentialOverlaps = findPotentialOverlaps(distribution, diameters)
% 查找可能的交叠峰区域
    % 初始化
    potentialOverlaps = [];
    
    % 平滑处理以减少噪声影响
    smoothDist = smoothdata(distribution, 'gaussian', 5);
    
    % 计算二阶导数
    d2 = diff(diff(smoothDist));
    d2 = [0, d2, 0];  % 补齐长度
    
    % 检测曲率变化区域（可能是交叠峰）
    for i = 5:length(d2)-5
        % 寻找曲率有显著变化的区域
        if abs(d2(i)) > 2*std(d2) && distribution(i) > max(distribution)*0.1
            % 检查是否有肩峰特征
            leftSlope = smoothDist(i) - smoothDist(max(i-4, 1));
            rightSlope = smoothDist(min(i+4, end)) - smoothDist(i);
            
            if (leftSlope * rightSlope < 0) || abs(leftSlope/rightSlope) > 2 || abs(rightSlope/leftSlope) > 2
                % 可能是交叠区
                potentialOverlaps = [potentialOverlaps, i];
            end
        end
    end
    
    % 合并太接近的区域
    if length(potentialOverlaps) > 1
        i = 1;
        while i < length(potentialOverlaps)
            if potentialOverlaps(i+1) - potentialOverlaps(i) <= 5
                potentialOverlaps(i) = round(mean([potentialOverlaps(i), potentialOverlaps(i+1)]));
                potentialOverlaps(i+1) = [];
            else
                i = i + 1;
            end
        end
    end
end

function overlapParams = adjustParamsForOverlap(initialParams, numPeaks, overlapIdx, diameters)
% 根据交叠区域调整参数
    overlapParams = initialParams;
    
    % 确定交叠区域的直径
    overlapDiameter = diameters(overlapIdx);
    logOverlapDiameter = log(overlapDiameter);
    
    % 找到距离交叠区域最近的峰
    closestPeaks = [];
    distances = [];
    
    for i = 1:numPeaks
        % 获取峰的中心位置
        peakCenter = exp(initialParams(3*i-1)); 
        distance = abs(peakCenter - overlapDiameter);
        closestPeaks = [closestPeaks, i];
        distances = [distances, distance];
    end
    
    % 按距离排序
    [~, idx] = sort(distances);
    closestPeaks = closestPeaks(idx);
    
    % 选取最近的两个峰进行调整
    for i = 1:min(2, length(closestPeaks))
        peakIdx = closestPeaks(i);
        
        % 修改峰的中心，使其更靠近交叠区域
        currentCenter = initialParams(3*peakIdx-1);
        newCenter = currentCenter + 0.2 * (logOverlapDiameter - currentCenter);
        overlapParams(3*peakIdx-1) = newCenter;
        
        % 减小峰的宽度，使交叠更清晰
        currentWidth = initialParams(3*peakIdx);
        overlapParams(3*peakIdx) = currentWidth * 0.85;
        
        % 增大振幅
        currentAmp = initialParams(3*peakIdx-2);
        overlapParams(3*peakIdx-2) = currentAmp * 1.15;
    end
end

function qualityScore = evaluateFitQuality(fitResult, gof, distribution, diameters, numPeaks)
% 评估拟合质量，特别关注交叠峰的分辨能力
    % 基础得分为RMSE
    qualityScore = gof.rmse;
    
    % 获取拟合曲线
    fittedVals = fitResult(diameters);
    
    % 计算残差
    residuals = distribution - fittedVals;
    
    % 惩罚过大或不合理的参数
    coeffs = coeffvalues(fitResult);
    for i = 1:numPeaks
        amp = coeffs(3*i-2);
        center = coeffs(3*i-1);
        width = coeffs(3*i);
        
        % 惩罚过宽的峰
        if width > 1.0
            qualityScore = qualityScore * (1 + 0.2 * (width - 1));
        end
        
        % 惩罚过小的峰
        if amp < max(distribution) * 0.05
            qualityScore = qualityScore * 1.2;
        end
        
        % 惩罚靠近边缘的峰
        if exp(center) < min(diameters)*1.2 || exp(center) > max(diameters)*0.8
            qualityScore = qualityScore * 1.3;
        end
    end
    
    % 检查残差模式 - 系统性残差表明拟合不佳
    residualCorr = xcorr(residuals, 5, 'coeff');
    if max(abs(residualCorr(6:end))) > 0.5
        qualityScore = qualityScore * 1.2;
    end
    
    % 对交叠峰区域的拟合质量额外评估
    if numPeaks >= 2
        % 检查每对相邻峰之间的区域
        coeffs = coeffvalues(fitResult);
        centers = zeros(numPeaks, 1);
        
        for i = 1:numPeaks
            centers(i) = exp(coeffs(3*i-1));
        end
        
        % 对每对相邻峰之间的区域进行检查
        [sortedCenters, idx] = sort(centers);
        for i = 1:numPeaks-1
            peak1 = idx(i);
            peak2 = idx(i+1);
            
            % 计算两峰之间的区域
            midPoint = sqrt(sortedCenters(i) * sortedCenters(i+1));
            [~, midIdx] = min(abs(diameters - midPoint));
            
            % 检查该区域的拟合质量
            localRange = max(1, midIdx-5):min(length(diameters), midIdx+5);
            localResiduals = residuals(localRange);
            localRMSE = sqrt(mean(localResiduals.^2));
            
            % 如果交叠区域拟合特别差，增加惩罚
            if localRMSE > gof.rmse * 1.5
                qualityScore = qualityScore * (1 + 0.3 * (localRMSE/gof.rmse - 1));
            end
        end
    end
end

function fitResult = postProcessFit(fitResult, distribution, diameters, numPeaks)
% 后处理拟合结果，修正明显不合理的参数
    if numPeaks <= 1
        return;  % 单峰不需要后处理
    end
    
    coeffs = coeffvalues(fitResult);
    
    % 检查峰位是否合理
    centers = zeros(numPeaks, 1);
    amps = zeros(numPeaks, 1);
    for i = 1:numPeaks
        centers(i) = exp(coeffs(3*i-1));
        amps(i) = coeffs(3*i-2);
    end
    
    % 检查峰是否在合理范围内
    validRangeLow = min(diameters) * 0.9;
    validRangeHigh = max(diameters) * 1.1;
    
    for i = 1:numPeaks
        if centers(i) < validRangeLow || centers(i) > validRangeHigh || amps(i) < 0
            % 此峰不合理，需要修正
            warning('发现不合理的峰位，进行修正');
            
            % 寻找分布中的最大值作为替代
            [maxVal, maxIdx] = max(distribution);
            newCenter = log(diameters(maxIdx));
            
            % 修改拟合结果
            newCoeffs = coeffs;
            newCoeffs(3*i-2) = maxVal * 0.8;  % 振幅
            newCoeffs(3*i-1) = newCenter;      % 对数中心
            
            % 尝试重新拟合
            try
                fitType = fittype(formula(fitResult), 'independent', 'x');
                options = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', newCoeffs);
                fitResult = fit(diameters', distribution', fitType, options);
            catch
                warning('后处理拟合失败，保持原拟合结果');
            end
            
            break; % 只修正一个参数，避免过度修改
        end
    end
end

function fitType = createOverlappingPeakFitType(numPeaks)
% 为交叠峰创建特殊的拟合类型
    % 这个函数创建一个修改版的拟合类型，其中两个峰可能共享参数
    if numPeaks == 2
        % 两峰交叠的特殊情况：允许共享宽度
        fitType = fittype('a1*exp(-((log(x)-b1)^2)/(2*c^2)) + a2*exp(-((log(x)-b2)^2)/(2*c^2))', ...
                         'independent', 'x', ...
                         'coefficients', {'a1', 'b1', 'a2', 'b2', 'c'});
    else
        % 默认回到标准多峰模型
        expression = '';
        for i = 1:numPeaks
            if i > 1
                expression = [expression ' + '];
            end
            expression = [expression sprintf('a%d*exp(-((log(x)-b%d)^2)/(2*c%d^2))', i, i, i)];
        end
        fitType = fittype(expression, 'independent', 'x');
    end
end

function params = estimateOverlappingParams(distribution, diameters)
% 为交叠峰拟合估计参数
    % 假设这是两峰交叠的情况
    
    % 首先尝试查找两个主要区域
    smoothDist = smoothdata(distribution, 'gaussian', 5);
    [pks, locs] = findpeaks(smoothDist, 'MinPeakProminence', max(smoothDist)*0.05);
    
    if length(pks) >= 2
        % 取两个最高峰
        [~, idx] = sort(pks, 'descend');
        locs = locs(idx(1:2));
        pks = pks(idx(1:2));
    else
        % 找不到两个峰，尝试寻找肩峰
        d2 = diff(diff(smoothDist));
        d2 = [0, d2, 0];
        [~, inflectionPoints] = findpeaks(abs(d2), 'MinPeakProminence', max(abs(d2))*0.1);
        
        if ~isempty(inflectionPoints) && ~isempty(locs)
            % 使用现有峰和曲率变化点
            secondLoc = inflectionPoints(1);
            secondPk = distribution(secondLoc) * 0.9;
            locs = [locs(1); secondLoc];
            pks = [pks(1); secondPk];
        else
            % 强制分为两部分
            [maxVal, maxIdx] = max(distribution);
            locs = [maxIdx; round(maxIdx * 0.7)];
            pks = [maxVal; maxVal * 0.5];
        end
    end
    
    % 确保有两个位置
    if length(locs) < 2
        locs = [locs; round(length(distribution)/2)];
        pks = [pks; max(distribution) * 0.5];
    end
    
    % 排序，确保位置升序
    [locs, idx] = sort(locs);
    pks = pks(idx);
    
    % 估计共享的宽度参数
    totalWidth = abs(log(diameters(locs(2))) - log(diameters(locs(1))));
    c = totalWidth / 2; 
    
    % 构建参数数组 [a1, b1, a2, b2, c]
    params = [pks(1), log(diameters(locs(1))), pks(2), log(diameters(locs(2))), c];
end

function peakParams = extractPeakParameters(fitResult, numPeaks)
% Extract peak parameters from the fit result
    peakParams = struct('amplitude', cell(numPeaks, 1), ...
                        'center', cell(numPeaks, 1), ...
                        'width', cell(numPeaks, 1));
                    
    coeffNames = coeffnames(fitResult);
    coeffValues = coeffvalues(fitResult);
    
    for i = 1:numPeaks
        % Extract a (amplitude), b (log center), c (width)
        aIdx = find(strcmp(coeffNames, sprintf('a%d', i)));
        bIdx = find(strcmp(coeffNames, sprintf('b%d', i)));
        cIdx = find(strcmp(coeffNames, sprintf('c%d', i)));
        
        peakParams(i).amplitude = coeffValues(aIdx);
        peakParams(i).center = exp(coeffValues(bIdx));  % Convert from log
        peakParams(i).width = coeffValues(cIdx);
    end
end

function analyzeParameterEvolution(fitResults, timePoints)
% Analyze how peak parameters evolve over time
    numTimePoints = length(fitResults);
    
    % Extract consistent peaks across time
    % This is challenging because the number of peaks may vary
    % For a simple approach, we'll track the most prominent peaks
    
    % Example visualization
    figure;
    
    % Plot evolution of peak centers
    subplot(3,1,1);
    hold on;
    for t = 1:numTimePoints
        for p = 1:length(fitResults(t).peaks)
            plot(timePoints(t), fitResults(t).peaks(p).center, 'ko');
        end
    end
    xlabel('Time');
    ylabel('Peak Center (nm)');
    title('Evolution of Peak Centers');
    set(gca, 'YScale', 'log');
    
    % Plot evolution of peak amplitudes
    subplot(3,1,2);
    hold on;
    for t = 1:numTimePoints
        for p = 1:length(fitResults(t).peaks)
            plot(timePoints(t), fitResults(t).peaks(p).amplitude, 'bo');
        end
    end
    xlabel('Time');
    ylabel('Peak Amplitude');
    title('Evolution of Peak Amplitudes');
    
    % Plot evolution of peak widths
    subplot(3,1,3);
    hold on;
    for t = 1:numTimePoints
        for p = 1:length(fitResults(t).peaks)
            plot(timePoints(t), fitResults(t).peaks(p).width, 'ro');
        end
    end
    xlabel('Time');
    ylabel('Peak Width');
    title('Evolution of Peak Widths');
end