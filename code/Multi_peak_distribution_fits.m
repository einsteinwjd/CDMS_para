function [fitResults, fittedDistributions] = Multi_peak_distribution_fits(simulatedPN, sim_sizebin)
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
fitResults = struct('time', cell(numTimePoints, 1), 'peaks', cell(numTimePoints, 1));
fittedDistributions = zeros(numTimePoints, numSizeBins);

% Fit options
options = fitoptions('Method', 'NonlinearLeastSquares',...
                    'Lower', [],...
                    'Upper', [],...
                    'StartPoint', [],...
                    'MaxIter', 1000,...
                    'TolFun', 1e-8,...
                    'TolX', 1e-8);

% Process each time point
for t = 1:numTimePoints
    % Current particle size distribution
    currentDist = particleCounts(t,:);
    
    % Determine number of peaks (could be automated with peak finding)
    numPeaks = estimateNumberOfPeaks(currentDist, diameters);
    
    % Create fitting function with appropriate number of peaks
    % Each peak modeled as lognormal: a*exp(-((log(x)-b)^2)/(2*c^2))
    fitType = createMultiPeakFitType(numPeaks);
    
    % Set initial parameters based on distribution
    initialParams = estimateInitialParameters(currentDist, diameters, numPeaks);
    options.StartPoint = initialParams;
    
    % Perform the fit
    [fitResult, gof] = fit(diameters', currentDist', fitType, options);
    
    % Store results
    fitResults(t).time = timePoints(t);
    fitResults(t).peaks = extractPeakParameters(fitResult, numPeaks);
    fitResults(t).gof = gof;
    
    % Store fitted distribution
    fittedDistributions(t,:) = fitResult(diameters);
end

% Analyze time evolution of parameters
analyzeParameterEvolution(fitResults, timePoints);

end

function numPeaks = estimateNumberOfPeaks(distribution, diameters)
% 使用更鲁棒的方法估计分布中的峰值数量
    % 使用不同的平滑宽度来捕捉不同尺度的特征
    smoothDist = smoothdata(distribution, 'gaussian', 5);
    
    % 使用对数横坐标找峰值，更适合粒径分布
    logDiameters = log10(diameters);
    
    % 使用相对高度和宽度标准
    maxHeight = max(smoothDist);
    minHeight = max(0.05 * maxHeight, median(smoothDist) * 0.5);
    
    % 找到所有可能的峰值
    [pks, locs, w, p] = findpeaks(smoothDist, 'MinPeakHeight', minHeight, ...
                                  'MinPeakProminence', maxHeight*0.08, ...
                                  'MinPeakDistance', 3);
    
    % 根据峰值特征过滤（宽度、突出度）
    validPeaks = (p > maxHeight*0.08) & (w > 1);
    pks = pks(validPeaks);
    locs = locs(validPeaks);
    
    numPeaks = length(pks);
    
    % 保证至少有一个峰值，最多4个
    numPeaks = max(1, min(numPeaks, 4));
    
    % 可视化峰值检测结果（调试时使用）
    % figure; plot(diameters, distribution, 'b-'); hold on;
    % plot(diameters, smoothDist, 'k-');
    % if ~isempty(locs)
    %     plot(diameters(locs), pks, 'ro', 'MarkerFaceColor', 'r');
    % end
    % set(gca, 'XScale', 'log');
    % title(['检测到 ', num2str(numPeaks), ' 个峰值']);
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
% Estimate initial parameters for fitting
    params = [];
    
    % Smooth the distribution for better peak finding
    smoothDist = smoothdata(distribution, 'gaussian', 5);
    
    % Find peaks
    [pks, locs] = findpeaks(smoothDist, 'NPeaks', numPeaks, 'SortStr', 'descend');
    
    % If peaks not found, create evenly spaced peaks
    if length(pks) < numPeaks
        pks = linspace(max(distribution)*0.8, max(distribution)*0.2, numPeaks);
        locs = round(linspace(1, length(distribution), numPeaks));
    end
    
    % For each peak, estimate a, b, c parameters
    for i = 1:numPeaks
        if i <= length(pks)
            a = pks(i);  % Height
            b = log(diameters(locs(i)));  % Log of center diameter
            
            % Estimate width by finding half-max width
            halfMax = a/2;
            leftIdx = find(smoothDist(1:locs(i)) <= halfMax, 1, 'last');
            rightIdx = locs(i) + find(smoothDist(locs(i):end) <= halfMax, 1, 'first') - 1;
            
            if isempty(leftIdx), leftIdx = 1; end
            if isempty(rightIdx), rightIdx = length(diameters); end
            
            width = log(diameters(rightIdx)) - log(diameters(leftIdx));
            c = width/2.355;  % Convert FWHM to standard deviation
        else
            % Fallback values
            a = max(distribution)*0.5;
            b = log(diameters(round(length(diameters)/2)));
            c = 0.5;
        end
        
        % Ensure c is positive and not too small
        c = max(c, 0.1);
        
        params = [params, a, b, c];
    end
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