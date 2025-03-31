function flux = Jx_cal(simulatedPN, sizeGrid, targetSize)
% Jx_cal - Calculate the flux of particles growing beyond a specified size
% 
% Input:
%   simulatedPN - Structure containing simulated particle number data
%   sizeGrid - Array of particle size bins (nm)
%   targetSize - Target size for flux calculation (nm)
%
% Output:
%   flux - Flux of particles growing beyond the target size (particles/cm³/s)

    % Extract necessary data
    time_point = simulatedPN.Time;  % Time in seconds

    N = simulatedPN{:,:};  % Particle number concentration
    
    % Find the bin index corresponding to particles >= targetSize
    [~, idxTarget] = min(abs(sizeGrid - targetSize));
    if sizeGrid(idxTarget) < targetSize && idxTarget < length(sizeGrid)
        idxTarget = idxTarget + 1;
    end
    
    % Calculate total particles above target size at each time
    N_aboveTarget = sum(N(:, idxTarget:end), 2);
    
    % Calculate flux (rate of change)
    flux = zeros(size(time_point));
    for i = 2:length(time_point)
        dt = seconds(time_point(i) - time_point(i-1));
        flux(i) = (N_aboveTarget(i) - N_aboveTarget(i-1)) / dt;
    end
    % Apply time smoothing using moving average
    windowSize = 5;  % Number of points in the moving average window
    flux = movmean(flux, windowSize);
    % Only consider positive flux (growth across boundary)
    flux(flux < 0) = 0;
end