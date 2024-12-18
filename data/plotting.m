close all; clc;clear

bigImg = im2gray(imread('itu_map_square.jpg'));
bigImg = im2double(bigImg); % for safe interpolation if needed

% Dimensions of the camera image
camWidth = 320;
camHeight = 240;

% UAV position and yaw angle
yaw = 60; % degrees
UAV_loc = [1500, 1500];
UAV_startx = round(UAV_loc(1) - camWidth/2);
UAV_starty = round(UAV_loc(2) - camHeight/2);
UAV_image = imrotate(bigImg(UAV_starty:UAV_starty+camHeight-1, UAV_startx:UAV_startx+camWidth-1),yaw);

% Detect Features of big image
key_points = detectORBFeatures(bigImg,"NumLevels",1);
n_downsample = 10;
figure
imshow(bigImg)
hold on
plot(key_points.Location(1:n_downsample:end,1),key_points.Location(1:n_downsample:end,2),'+g')

% Add text to the bottom-right corner
xlim_values = xlim; % Get the x-axis limits
ylim_values = ylim; % Get the y-axis limits

x_position1 = xlim_values(2) - 0.05 * (xlim_values(2) - xlim_values(1)); % Slight padding from right
y_position1 = ylim_values(2) - 0.05 * (ylim_values(2) - ylim_values(1)); % Slight padding from bottom

y_position2 = ylim_values(2) - 0.02 * (ylim_values(2) - ylim_values(1)); % Slight padding from bottom


text(x_position1, y_position1, ['Detected Features: ', num2str(key_points.Count) ], 'Color', 'red', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 12,'BackgroundColor','Black');

text(x_position1, y_position2, ['Plotted Features: ', num2str(round(key_points.Count/n_downsample)) ], 'Color', 'red', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 12,'BackgroundColor','Black');

title('Feature Detection using ORB Algorithm(Hand-made Local Feature Detector)')

%%
bigImg = (imread('itu_map_square.jpg'));

imshow(bigImg)

title('High Resolution Google Earth Satellite Image')
