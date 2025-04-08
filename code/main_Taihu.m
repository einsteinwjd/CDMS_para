clear; close all;

% path define
path_define;

% exe
obs_to_timetable;
modeldata_to_timetable;
% timetable_plot; % plot size distribution for specific time
for target_day = 1:6
    Moments_fits(target_day,'path_define'); % Moments method to fit diamter
    Moments_fits_obs(target_day,'path_define');

end
% multi_peak_fit; % multi-peak fit method, using Multi_peak_distribution_fits functions.