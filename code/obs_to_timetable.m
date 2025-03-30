clear; close all;
%% Load Taihu NPF data
% This part loads the Taihu NPF data from the MAT file and displays the
% Read the MAT file containing Taihu NPF data
path_define;

F0_path = [F0_folder,'Taihu_NPFData.mat']
% Check if file exists
if ~exist(F0_path, 'file')
    error('Data file not found: %s', F0_path);
end

% Load the MAT file
try
    data = load(F0_path);
catch ME
    error('Error loading data file: %s', ME.message);
end

% Display the contents of the loaded data
disp('Available variables in the loaded data:');
disp(fieldnames(data));

%% Extract the data
% generate time labels. The data is sampled every 12 minutes
% the date is extracted from data.date. The format is 'yyyymmdd'
% Convert date strings to datetime objects
dates = datetime(string(data.date), 'InputFormat', 'yyyyMMdd');

% Create time vector for one day (from 0 to 24 hours, 12-minute intervals)
minutes_per_day = 24 * 60;
time_intervals = 12; % 12 minutes
time_points = minutes_per_day / time_intervals;
time_vector = (0:time_points-1) * (time_intervals/60); % in hours

% Create a complete datetime array for all measurements
[D, T] = meshgrid(dates, time_vector);
time_labels = D + hours(T);
time_labels = time_labels(:); % Convert to column vector

% Extract data to timetable
% Concatenate data from cells
sa_combined = vertcat(data.sa{:});
dma_combined = vertcat(data.dma{:});
org_combined = vertcat(data.org{:});

% Create variable names for org columns
org_names = arrayfun(@(x) ['vbs' num2str(x)], 1:9, 'UniformOutput', false);

% Create timetable
tt_precursors = timetable(time_labels, sa_combined, dma_combined, org_combined(:,1), org_combined(:,2), org_combined(:,3), org_combined(:,4), org_combined(:,5), org_combined(:,6), org_combined(:,7), org_combined(:,8), org_combined(:,9));

% Set variable names
tt_precursors.Properties.VariableNames = {'sa', 'dma', org_names{:}};

% generate a new timetable for data.psd
% the time is time_labels
% the psd data is in data.psd
% the psd data is in a cell array. Each cell contains a matrix with 159 columns
% the column is the diameter in nm, with size determined by data.size_vec
% the row is the time
% Concatenate PSD data from cells
psd_combined = vertcat(data.psd{:});

% Create variable names for PSD columns using size_vec
psd_names = arrayfun(@(x) ['d' num2str(x, '%.1f')], data.size_vec, 'UniformOutput', false);

% Create timetable for PSD data
% Convert time and PSD data into timetable more efficiently
tt_psd = array2timetable(psd_combined, 'RowTimes', time_labels);

% Set variable names
tt_psd.Properties.VariableNames = psd_names;

%% save data
% save the tt_psd and tt_precurosrs to a MAT file in F1 folder
% Create F1 directory if it doesn't exist
if ~exist(F1_folder, 'dir')
    mkdir(F1_folder);
end

% Save the timetables
save(fullfile(F1_folder, 'timetable_data.mat'), 'tt_psd', 'tt_precursors');
