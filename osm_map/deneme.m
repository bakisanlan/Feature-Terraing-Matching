name = "openstreetmap";
url = "https://basemaps.arcgis.com/arcgis/rest/services/OpenStreetMap_v2/VectorTileServer/tile/${z}/${y}/${x}.pbf";
style = "https://www.arcgis.com/sharing/rest/content/items/3e1a00aeae81496587988075fe529f71/resources/styles/root.json?f=pjson";
addCustomBasemap(name,url,Style=style)

%%

latlim = [41.09831 41.10966];
lonlim = [29.01378 29.03970];
% [A,RA,attrib] = readBasemapImage(name,latlim,lonlim);

[A,RA,attribA] = readBasemapImage("satellite",latlim,lonlim);

%%

figure
mapshow(A,RA)
axis off

plt = gca;
exportgraphics(plt,'itu_sat.jpg','Resolution',600)


%% crop and rotate image
I = imread('itu_sat.jpg');
[J,rect] = imcrop(I);
J = imrotate(J,-45,'bilinear','crop');

imwrite(J,'croped_turned_itu.jpg')

%% detect feature 
tic
%org image feature detection through SURF
RGB1 = imread('itu_sat.jpg');
I1 = rgb2gray(RGB1);
points1 = detectSURFFeatures(I1);
[features1,valid_points1] = extractFeatures(I1,points1);
toc

%cropped and turned image feature detection through SURF
tic
RGB2 = imread('croped_turned_itu.jpg');
I2 = rgb2gray(RGB2);
points2 = detectSURFFeatures(I2);
[features2,valid_points2] = extractFeatures(I2,points2);
toc

%% match feature
tic 
indexPairs = matchFeatures(features1,features2);

matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);

% showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2);
toc
%% exclude outliers using RANSAC

tic
[tform,inlierIdx] = estgeotform2d(matchedPoints2,matchedPoints1,"similarity");
inlierPtsDistorted = matchedPoints2(inlierIdx,:);
inlierPtsOriginal  = matchedPoints1(inlierIdx,:);
 
% figure
% showMatchedFeatures(I1,I2,inlierPtsOriginal,inlierPtsDistorted)
title("Matched Inlier Points")
toc

%% recover distorted image
outputView = imref2d(size(I1));
Ir = imwarp(I2,tform,"OutputView",outputView);
 
imshow(Ir) 
title("Recovered Image")
hold on
imshow(I1)