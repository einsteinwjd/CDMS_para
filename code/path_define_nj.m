% Define folder paths
F0_folder = '../data/F0_nj/';  % Path to F0 folder
F1_folder = '../data/F1_nj/';  % Path to F1 folder
F2_folder = '../data/F2_nj/';  % Path to F1 folder
Fig_folder = '../figs_nj/';  % Path to figures folder

% Verify if folders exist
if ~exist(F0_folder, 'dir')
    warning('F0 folder does not exist at: %s', F0_folder);
end

if ~exist(F1_folder, 'dir')
    warning('F1 folder does not exist at: %s', F1_folder);
end