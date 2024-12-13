%% detect feature
tic
%org image feature detection through SURF
RGB1 = imread('croped_arc.jpg');
I1 = rgb2gray(RGB1);
points1 = detectSURFFeatures(I1);
[features1,valid_points1] = extractFeatures(I1,points1);
toc

%cropped and turned image feature detection through SURF

for i=0:100

    % tic
    RGB2 = imread(['crop_imgs\croped_',num2str(i),'.jpg']);
    I2 = rgb2gray(RGB2);
    points2 = detectSURFFeatures(I2);
    [features2,valid_points2] = extractFeatures(I2,points2);
    % toc

    %% match feature
    % tic
    indexPairs = matchFeatures(features1,features2);

    matchedPoints1 = valid_points1(indexPairs(:,1),:);
    matchedPoints2 = valid_points2(indexPairs(:,2),:);

    if size(matchedPoints1.Location,1) < 2
        continue
    else
        %% exclude outliers using RANSAC

        % tic
        [tform,inlierIdx,stat] = estgeotform2d(matchedPoints2,matchedPoints1,"similarity");

        if stat ~= 0
            continue
        else
            inlierPtsDistorted = matchedPoints2(inlierIdx,:);
            inlierPtsOriginal  = matchedPoints1(inlierIdx,:);


            % figure
            figure
            showMatchedFeatures(I1,I2,inlierPtsOriginal,inlierPtsDistorted)
            title("Matched Inlier Points")
            % toc
        end
    end

end
