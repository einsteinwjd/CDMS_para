clear; close all;

% path define
path_define_nj;

% exe
obs_to_timetable_nj;
modeldata_to_timetable_nj;
% timetable_plot_nj;
for target_day = 1:3
    Moments_fits(target_day,'path_define_nj'); % Moments method to fit diamter
    Moments_fits_obs(target_day,'path_define_nj');

end
% multi_peak_fit; % multi-peak fit method, using Multi_peak_distribution_fits functions.