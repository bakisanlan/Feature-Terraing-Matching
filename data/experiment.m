% Load the big aerial image
close all; clc;clear
rng(5);
bigImg = im2gray(imread('itu_map_square.jpg'));
bigImg = im2double(bigImg); % for safe interpolation if needed


N = 500;

% Dimensions of the camera image
camWidth = 320;
camHeight = 240;

% UAV position and yaw angle
l = 1000;
u = 2000;
Center = l + (u-l)*rand(N,2);
partPosXCell = num2cell(Center(:,1));
partPosYCell = num2cell(Center(:,2));


% xCenter = 2500;
% yCenter = 2000;
yaw = 20; % degrees
UAV_loc = [1500, 1500];
UAV_startx = round(UAV_loc(1) - camWidth/2);
UAV_starty = round(UAV_loc(2) - camHeight/2);

UAV_image = imrotate(bigImg(UAV_starty:UAV_starty+camHeight-1, UAV_startx:UAV_startx+camWidth-1),yaw,'crop');

% Convert yaw to radians
theta = deg2rad(yaw);

% Prepare coordinate grids for the UAV camera frame
[uGrid, vGrid] = meshgrid(1:camWidth, 1:camHeight);

% Shift origin to center of the UAV image
uCenter = (camWidth+1)/2;
vCenter = (camHeight+1)/2;
uPrime = uGrid - uCenter;
vPrime = vGrid - vCenter;

% Apply rotation by yaw angle
% Rotation matrix: [cos(theta) -sin(theta); sin(theta) cos(theta)]
xRel =  uPrime*cos(theta) - vPrime*sin(theta);
yRel =  uPrime*sin(theta) + vPrime*cos(theta);

% Translate to big image coordinates
% xImg = xCenter + xRel;
% yImg = yCenter + yRel;

tic
xImgCell = cellfun(@(x) xRel+x, partPosXCell,'UniformOutput',false);
yImgCell = cellfun(@(y) yRel+y, partPosYCell,'UniformOutput',false);
toc

% Make sure coordinates are within image bounds (simple clamping)
[H, W, C] = size(bigImg);
% xImg(xImg < 1) = 1; xImg(xImg > W) = W;
% yImg(yImg < 1) = 1; yImg(yImg > H) = H;

% Interpolate from bigImg
% If RGB: separate channels
if C == 3
    part_image = zeros(camHeight, camWidth, 3);
    for c = 1:3
        % part_image(:,:,c) = interp2(1:W, 1:H, bigImg(:,:,c), xImg, yImg, 'linear');

        part_image = cellfun(@(x,y) interp2(1:W, 1:H, bigImg(:,:,c), x, y, 'linear'),xImgCell,yImgCell,'UniformOutput',false);
    end
else
    % Grayscale
    % part_image = interp2(1:W, 1:H, bigImg, xImg, yImg, 'linear');
    part_image = cellfun(@(x,y) interp2(1:W, 1:H, bigImg, x, y, 'linear'),xImgCell,yImgCell,'UniformOutput',false);

end
toc
% Display the simulated UAV image
% figure
% imshow(part_image);
% gray1 = im2gray(part_image);
% title(sprintf('Simulated UAV Image at (%.0f,%.0f), Yaw=%.0f°', xCenter, yCenter, yaw));

%%
close all; clc;clear
% rng(5);
bigImg = im2gray(imread('itu_map_square.jpg'));
bigImg = im2double(bigImg); % for safe interpolation if needed

N = 10;

% Dimensions of the camera image
camWidth = 320;
camHeight = 240;

% UAV position and yaw angle
l = 1000;
u = 2000;
Center = l + (u-l)*rand(N,2);
xStartCell = num2cell(round(Center(:,1) - camWidth/2));
yStartCell = num2cell(round(Center(:,2) - camHeight/2));

yaw = 0; % degrees
UAV_loc = [1500, 1500];
UAV_startx = round(UAV_loc(1) - camWidth/2);
UAV_starty = round(UAV_loc(2) - camHeight/2);
UAV_image = imrotate(bigImg(UAV_starty:UAV_starty+camHeight-1, UAV_startx:UAV_startx+camWidth-1),yaw);

% tic
part_image = cellfun(@(x,y) imrotate(bigImg(y:y+camHeight-1, x:x+camWidth-1, :), yaw),xStartCell,yStartCell,'UniformOutput',false);
% toc

%feature extraction
tic
[features1,valid_points1] = extractFeatures(UAV_image,detectORBFeatures(UAV_image,"NumLevels",1));
toc

tic
[features2,valid_points2] = cellfun(@(x) extractFeatures(x,detectORBFeatures(x,"NumLevels",1)),part_image,'UniformOutput',false);
toc

tic
indexPairs = cellfun(@(x) matchFeatures(features1,x,"Method","Exhaustive"),features2,'UniformOutput',false);
toc

tic
[~,inlierIdx,~] = cellfun(@(PartValid,idxPair) estgeotform2d(PartValid(idxPair(:,2),:),valid_points1(idxPair(:,1),:),"similarity"),valid_points2,indexPairs,'UniformOutput',false);
toc

tic
numMatchedFeature = cellfun(@(x) sum(x),inlierIdx);
toc

%

% tic
% for i=1:N
%     [features2,valid_points2] = extractFeatures(part_image{i},detectORBFeatures(part_image{i},"NumLevels",1));
% 
%     indexPairs = matchFeatures(features1,features2);
% 
%     [~,inlierIdx,~] = estgeotform2d(valid_points2(indexPairs(:,2),:),valid_points1(indexPairs(:,1),:),"similarity");
% end
% toc


%%


% Show the large base map
figure;
imshow(imread('itu_map_square.jpg'));
hold on;

for i = 1:N
    xCenter = Center(i,1);
    yCenter = Center(i,2);

    % Compute the top-left corner of the image so that its center is at (xCenter, yCenter).
    xStart = round(xCenter - camWidth/2);
    yStart = round(yCenter - camHeight/2);

    % % Ensure we don’t exceed map boundaries (optional):
    % [H, W, ~] = size(bigImg);
    % if xStart < 1, xStart = 1; end
    % if yStart < 1, yStart = 1; end
    % if xStart + camWidth - 1 > W, xStart = W - camWidth + 1; end
    % if yStart + camHeight - 1 > H, yStart = H - camHeight + 1; end

    % Plot the image at the given location
    % Note: 'YData' is given top-to-bottom, so we specify [yStart yStart+camHeight-1].
    % Similarly, XData is left-to-right.
    image('XData', [xStart xStart+camWidth-1], 'YData', [yStart yStart+camHeight-1], 'CData', part_image{i}*255);

    matchedPoints1 = valid_points1(indexPairs{i}(:,1),:).Location;
    matchedPoints2 = valid_points2{i}(indexPairs{i}(:,2),:).Location;

    plot(matchedPoints2(:,1)+xStart,matchedPoints2(:,2)+yStart,'+r');
end

hold off;
title('All 500 images overlaid at their respective locations');



%%
xCenter = 1500;
yCenter = 2500;
camWidth = 640;
camHeight = 480;
xstart = xCenter - camWidth/2;
ystart = yCenter - camHeight/2;

bigImg = imread('itu_map_square.jpg');
subimage = bigImg(ystart:ystart+camHeight-1,xstart:xstart+camWidth-1,:);
yaw = 45;
% distortedimage = imresize(subimage,0.5); 
% distortedimage    = imrotate(subimage,yaw);
% gray2 = im2gray(distortedimage);
figure
imshow(subimage)

%%

figure
% points1 = detectSIFTFeatures(UAV_image,'NumLayersInOctave',4,'ContrastThreshold',0.0066667);
points1 = detectORBFeatures(UAV_image);
[features1,valid_points1] = extractFeatures(UAV_image,points1);

% points2 = detectSIFTFeatures(distortedimage,'NumLayersInOctave',4,'ContrastThreshold',0.0066667);
points2 = detectORBFeatures(distortedimage);

[features2,valid_points2] = extractFeatures(distortedimage,points2);

indexPairs = matchFeatures(features1,features2);

matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);

[tform,inlierIdx] = estgeotform2d(matchedPoints1,matchedPoints2,"similarity");

inlierPts1 = matchedPoints1(inlierIdx,:);
inlierPts2 = matchedPoints2(inlierIdx,:);
 
showMatchedFeatures(UAV_image,distortedimage,inlierPts1,inlierPts2,'montage',PlotOptions={'go','go','g-'})

xlim_values = xlim; % Get the x-axis limits
ylim_values = ylim; % Get the y-axis limits

x_position1 = xlim_values(1) + 3; % Slight padding from right
y_position1 = ylim_values(1) + 3; % Slight padding from bottom

y_position2 = ylim_values(2) - 0.02 * (ylim_values(2) - ylim_values(1)); % Slight padding from bottom


text(x_position1, y_position1,['Inlier Matched Features: ', num2str(sum(inlierIdx)),'/',num2str(length(inlierIdx))], 'Color', 'red', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 20,'BackgroundColor','Black');

