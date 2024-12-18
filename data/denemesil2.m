%% Setup
% Assume we have a large satellite image: bigImage (4000x4000)
% For demonstration, let's create a dummy image if you don't have one:
% bigImage = imread('YourSatelliteImage.png'); % real use
bigImage = (imread('itu_map_square.jpg'));

% UAV camera properties
cameraWidth = 320;
cameraHeight = 240;

% UAV position in the big image (pixels)
uavCenterX = 2000; 
uavCenterY = 2000;

% UAV orientation angles (in degrees)
yawAngle = 0;   % Rotation around the vertical (Z) axis
pitchAngle = 30;  % Small tilt forward/back (rotation about the Y-axis)
rollAngle = 30;   % Small tilt left/right (rotation about the X-axis)

% Scaling factor to simulate altitude difference
scaleFactor = 1.5; 

%% Step 1: Initial Crop
% We'll crop a larger region first, then rotate, project, and finally extract the 320x240 patch.

overshootFactor = 2; 
initialCropWidth = round(cameraWidth * scaleFactor * overshootFactor);
initialCropHeight = round(cameraHeight * scaleFactor * overshootFactor);

% Ensure crop is within image boundaries
startX = max(uavCenterX - initialCropWidth/2, 1);
startY = max(uavCenterY - initialCropHeight/2, 1);
endX = min(startX + initialCropWidth - 1, size(bigImage,2));
endY = min(startY + initialCropHeight - 1, size(bigImage,1));
startX = endX - initialCropWidth + 1; 
startY = endY - initialCropHeight + 1;

croppedImage = bigImage(startY:endY, startX:endX, :);

%% Step 2: Yaw Rotation
% Rotate about Z-axis (vertical) using imrotate
rotatedImage = imrotate(croppedImage, yawAngle, 'bicubic', 'crop');

%% Step 3: Projective Transform for Pitch & Roll
% Pitch and roll cause perspective distortion. We can model this as a projective transform.
% For small angles, we can approximate a homography that tilts the image.

% Convert angles to radians
pitchRad = deg2rad(pitchAngle);
rollRad = deg2rad(rollAngle);

% We will construct a homography that simulates these tilts.
% One simple approach is to consider how a flat scene would look if we rotate it in 3D.
% For a camera looking down (Z axis), a small pitch and roll introduce perspective foreshortening.

% A simplistic approach: 
% Let's assume the camera image plane is initially aligned with the ground.
% Pitch tilts along the camera's x-axis (left-right in image), causing vertical perspective distortion.
% Roll tilts along the camera's y-axis (up-down in image), causing horizontal perspective distortion.

% We can create a projective transform matrix H:
% Start with identity
H = eye(3);

% To approximate pitch and roll:
% - Roll around x-axis: will cause a vertical "shear" perspective.
% - Pitch around y-axis: will cause a horizontal "shear" perspective.

% A more physically motivated approach:
% If the camera tilts forward (pitch), distant objects appear higher/lower:
% We can simulate this by adding a small perspective term in the transform.
%
% Similarly, roll tilts the image side-to-side.

% Let's define a small factor related to angle. For small angles, 
% perspective distortions ~ tan(angle).
pxFactor = tan(pitchRad)*0.001; % Adjust these scaling factors as needed for realism
rlFactor = tan(rollRad)*0.001;

% Construct a projective transform that introduces perspective distortion.
% We'll do two steps: one for roll and one for pitch.
% For pitch: tilt in vertical direction (y-direction)
% Hpitch = [1     0      0;
%           0     1      pxFactor;
%           0     pxFactor 1];
% 
% % For roll: tilt in horizontal direction (x-direction)
% Hroll = [1    rlFactor 0;
%          rlFactor 1     0;
%          0       0     1];

theta = 10;
tm = [cosd(theta) -sind(theta) 0.001; ...
      sind(theta) cosd(theta) 0.01; ...
      0 0 1];
tform = projective2d(tm);

% Combine them: the final homography H = Hpitch * Hroll
% H = Hpitch * Hroll;

% Create a projective2d object
% tform = projective2d(H);

% Apply the projective transform
projectedImage = imwarp(rotatedImage, tform);

%% Step 4: Rescale the image (simulate altitude)
rescaledImage = imresize(projectedImage, 1/scaleFactor, 'bicubic');

%% Step 5: Final crop to camera size (320x240)
finalWidth = cameraWidth;
finalHeight = cameraHeight;

centerX = size(rescaledImage,2)/2;
centerY = size(rescaledImage,1)/2;

cropStartX = round(centerX - finalWidth/2);
cropStartY = round(centerY - finalHeight/2);

cropStartX = max(cropStartX, 1);
cropStartY = max(cropStartY, 1);

endX = min(cropStartX + finalWidth - 1, size(rescaledImage,2));
endY = min(cropStartY + finalHeight - 1, size(rescaledImage,1));

uavView = rescaledImage(cropStartY:endY, cropStartX:endX, :);

%% Step 6: Adjust lighting conditions
gammaValue = 1.2; 
adjustedImage = imadjust(uavView, [], [], gammaValue);

% Slight brightness increase
adjustedImage = adjustedImage + 0.05; 
adjustedImage(adjustedImage>1) = 1;

% You can also add noise or color adjustments here if desired.

%% Display results
figure;
subplot(1,3,1); imshow(bigImage); title('Original Satellite Image');
subplot(1,3,2); imshow(uavView); title('UAV View With Tilt (Before Lighting)');
subplot(1,3,3); imshow(adjustedImage); title('UAV View With Tilt (After Lighting)');
