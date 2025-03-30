% Define folder paths
F0_folder = '../data/F0/';  % Path to F0 folder
F1_folder = '../data/F1/';  % Path to F1 folder
F2_folder = '../data/F2/';  % Path to F1 folder


% Verify if folders exist
if ~exist(F0_folder, 'dir')
    warning('F0 folder does not exist at: %s', F0_folder);
end

if ~exist(F1_folder, 'dir')
    warning('F1 folder does not exist at: %s', F1_folder);
end