clear; close all;
%% load F1 folder timetable_data.mat
path_define;
load([F1_folder,'timetable_data.mat']);
load([F1_folder,'modeldata_to_timetable.mat']);
%% plot timetable

% Choose a specific time (e.g., first timestamp)
timepoint = 2320;

% Extract data for the chosen time
sizeData = simulatedPN{timepoint,:};

% Create figure with specified size
figure('Position', [100, 100, 800, 600], 'Color', 'white');

% Plot data with improved styling
p = plot(sim_sizebin, sizeData, 'LineWidth', 2, 'Color', [0.2, 0.4, 0.8]);

% Customize the axes
ax = gca;
set(ax, 'XScale', 'log', 'FontName', 'Arial', 'FontSize', 12, 'Box', 'on', 'LineWidth', 1.5);
set(ax, 'XGrid', 'on', 'YGrid', 'on', 'GridLineStyle', ':');
xlim([min(sim_sizebin) max(sim_sizebin)]);

% Add labels and title with improved formatting
xlabel('Particle Diameter (nm)', 'FontWeight', 'bold', 'FontSize', 14);
ylabel('Particle Number Concentration (#/cm³)', 'FontWeight', 'bold', 'FontSize', 14);
title(['Particle Size Distribution at ', datestr(timepoint)], 'FontWeight', 'bold', 'FontSize', 16);

% Add minor grid lines
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridAlpha = 0.25;

% Add a box around the plot
box on;

% Optionally add legend if needed
% legend('Measured Data', 'Location', 'best');

% Set tight figure margins
set(gcf, 'PaperPositionMode', 'auto');