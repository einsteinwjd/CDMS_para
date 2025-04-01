clear; close all;
path_define;

%% load model data from F0
files = dir(fullfile(F0_folder, 'output*.mat'));
data = cell(length(files), 1);

for i = 1:length(files)
    data{i} = load(fullfile(F0_folder, files(i).name));
end

%% Extract the data
% Concatenate outputPN from all files
simulatedPN = [];
for i = 1:length(data)
    simulatedPN = [simulatedPN; data{i}.outputPN];
end
sim_sizebin = data{1}.MP.dpSec; %nm

%% transform to timetable
% generate time vector rom 6:00 to 16:00 with 1 minute interval
% date is 20230723, 20230729, 20230730, 20230731, 20230803, 20230804

dates = [datetime(2023, 7, 23), datetime(2023, 7, 29), datetime(2023, 7, 30), datetime(2023, 7, 31), datetime(2023, 8, 3), datetime(2023, 8, 4)];
simulatedTime = [];
for i = 1:length(dates)
    simulatedTime = [simulatedTime, dates(i)+timeofday(datetime(0,0,0,8,0,0)):seconds(30):dates(i)+timeofday(datetime(0,0,0,18,0,0))];
end

% simulatedPN = simulatedPN';
% simulatedPN = simulatedPN(1:length(simulatedTime),:);
simulatedTime = simulatedTime(:);  % Convert to column vector
simulatedPN = array2timetable(simulatedPN, 'RowTimes', simulatedTime);

%% save to F1_folder
save([F1_folder, 'modeldata_to_timetable.mat'], 'simulatedPN','sim_sizebin');
disp('model data saved to F1_folder');