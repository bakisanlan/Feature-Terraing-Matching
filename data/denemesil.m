%% Setup
% Assume we have a large satellite image: bigImage (4000x4000)
% For demonstration, let's create a dummy image if you don't have one:
% bigImage = imread('YourSatelliteImage.png'); % real use
bigImage = im2gray(imread('itu_map_square.jpg'));


% UAV camera properties
cameraWidth = 320;
cameraHeight = 240;

% UAV position in the coordinate system of the big image
% Assuming you know the UAV’s "pixel" position in the big image:
uavCenterX = 2000; % example center X coordinate on big image
uavCenterY = 2000; % example center Y coordinate on big image

% UAV yaw angle (in degrees, clockwise)
yawAngle = 30; % for example

% Scaling factor to simulate altitude difference
% Suppose the satellite image represents ground at 1 meter/pixel 
% and UAV altitude + camera lens properties give a coverage 
% that corresponds to scaling. For example:
scaleFactor = 1.5; % If >1, we zoom out (to represent higher altitude), <1 zoom in

%% Step 1: Crop the big image at UAV position
% We want a region that, after rotation and scaling, results in a 320x240 image.
% To be safe, first we will extract a larger square patch from bigImage,
% rotate and scale it, then crop to desired 320x240.

% Determine the size of the patch to extract before rotation/rescaling.
% If final desired is 320x240, and we are scaling by factor and rotating,
% let's overshoot the initial crop to avoid edge cutoff. For instance:
overshootFactor = 2; 
initialCropWidth = round(cameraWidth * scaleFactor * overshootFactor);
initialCropHeight = round(cameraHeight * scaleFactor * overshootFactor);

% Ensure the crop does not go out of boundaries:
startX = max(uavCenterX - initialCropWidth/2, 1);
startY = max(uavCenterY - initialCropHeight/2, 1);

% Make sure the crop fits within the big image:
endX = min(startX + initialCropWidth - 1, size(bigImage,2));
endY = min(startY + initialCropHeight - 1, size(bigImage,1));

% Adjust startX, startY if endX, endY trimmed the rectangle:
startX = endX - initialCropWidth + 1;
startY = endY - initialCropHeight + 1;

croppedImage = bigImage(startY:endY, startX:endX, :);

%% Step 2: Rotate the image by the UAV's yaw angle
% Use imrotate (by default it pads with black, you can change 'crop' or interpolation method)
rotatedImage = imrotate(croppedImage, yawAngle, 'bicubic', 'crop');

%% Step 3: Rescale the image
% Rescale to simulate altitude difference. If scaleFactor > 1, you get a smaller view.
% If scaleFactor < 1, you get a more zoomed-in view.
rescaledImage = imresize(rotatedImage, 1/scaleFactor, 'bicubic');

%% Step 4: Final crop to the exact camera size (320x240)
% After rotation and scaling, the image might not match desired size exactly.
% We will crop it to the desired camera size. If the rescaled image is smaller,
% you might need a larger initial crop or smaller scaleFactor.
finalWidth = cameraWidth;
finalHeight = cameraHeight;

% Compute center of the rescaled image
centerX = size(rescaledImage,2)/2;
centerY = size(rescaledImage,1)/2;

cropStartX = round(centerX - finalWidth/2);
cropStartY = round(centerY - finalHeight/2);

% Ensure indexes are within image bounds
cropStartX = max(cropStartX, 1);
cropStartY = max(cropStartY, 1);

endX = min(cropStartX + finalWidth - 1, size(rescaledImage,2));
endY = min(cropStartY + finalHeight - 1, size(rescaledImage,1));

uavView = rescaledImage(cropStartY:endY, cropStartX:endX, :);

%% Step 5: Adjust lighting conditions
% The satellite image might look different from a UAV’s camera feed. 
% You might want to simulate different illumination, contrast, or add slight noise.

% Example adjustments:
% - Adjust gamma to simulate different lighting conditions
gammaValue = 1.2; % >1: brighter mid-tones, <1: darker mid-tones
adjustedImage = imadjust(uavView, [], [], gammaValue);

% - Add some random noise (to simulate camera sensor noise)
adjustedImage = imnoise(adjustedImage, 'gaussian', 0, 0.0005);

% - Adjust overall brightness or contrast
% For instance, increase brightness:
% adjustedImage = adjustedImage + 0.05; 
% adjustedImage(adjustedImage>1) = 1; % Ensure we don't exceed max

% - Adjust color slightly (simulate different white balance)
% Convert to HSV and shift hue/saturation if you wish:
% hsvImg = rgb2hsv(adjustedImage);
% hsvImg(:,:,2) = hsvImg(:,:,2) * 1.1; % Increase saturation slightly
% adjustedImage = hsv2rgb(hsvImg);

%% Step 6: Blur the image
ww = 3;
kernel = ones(ww)/ww^2;
adjustedImage = imfilter(adjustedImage,kernel);
%% Display results
figure;
subplot(1,3,1); imshow(croppedImage); title('Original Large Satellite Image');
subplot(1,3,2); imshow(uavView); title('UAV View Before Lighting Adjustment');
subplot(1,3,3); imshow(adjustedImage); title('UAV View After Lighting Adjustment');
