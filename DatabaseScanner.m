classdef DatabaseScanner < handle
    %DATABASESCANNER Class wil be used by State Estimators(PF) for
    %getting measurements of particles through scannig the offline
    %database map(Satellite Image).

   
    properties
        snapDim    %Snapped image dimension [W,H]
        AIM        %AerialImageModel object holder
        num_level  % ORB num level parameter

    end

    methods

        function obj = DatabaseScanner(varargin)
            %DatabaseScanner Construct an instance of this class

            % Call parent class constructor
            obj.snapDim = [320,240];
            obj.num_level = 1;

        end

        function partlikelihood = find_likelihood(obj,partImgCenterWorldPos,yaw,UAVImage)
            % find_likelihood function will be called by State Estimators(PF)
            % to get likelihood vaulues of particles. Likelihood values 
            % are calculated through image mathcing prosedure. 

            % Particles Image Cell Array that created by snapPartImage
            % function with using 2D position and yaw information of
            % particles
            PartImage = obj.snapPartImage(partImgCenterWorldPos,yaw);

            % Getting inlierIdx boolean(0 or 1) cell array
            inlierIdx = obj.ImageMatching(UAVImage,PartImage);

            % Counting how many inliers matched feature points have each
            % particles with using sum(x) function
            numMatchedFeature = cellfun(@(x) sum(x),inlierIdx);

            % For now we find likelihood = number of mathced inlier feature
            %points. Later we can use an exponential function to handle
            %beter results.
            partlikelihood = numMatchedFeature;
        end

        function UAVImage = snapUAVImage(obj,UAVWorldPos,yaw)
        
            % Rescale y component due to direction of y axis of image frame
            % is opposite with world frame
            UAVWorldPos(2) = -UAVWorldPos(2);

            % Finding UAV corresponding pixel location using meter/pixel
            % ratio
            UAVImgCenterpxPos = round(UAVWorldPos * (1/obj.AIM.mp));

            % Dimensions of the camera image
            w = obj.snapDim(1);
            h = obj.snapDim(2);

            % Cropping and rotating(might be unnecessary) UAV image
            xStart = UAVImgCenterpxPos(1) - w/2;
            yStart = UAVImgCenterpxPos(2) - h/2;
            UAVImage = imrotate(obj.AIM.I(yStart:yStart+h-1, xStart:xStart+w-1, :), yaw, 'crop');
    
        end

        function PartImage = snapPartImage(obj,partWorldPos,yaw)

            % Rescale y component due to direction of y axis of image frame
            partWorldPos(:,2) = -partWorldPos(:,2);

            % Convert (N,E)meters dimension to (u,v)px
            partImgCenterpxPos = round(partWorldPos * (1/obj.AIM.mp));

            % Dimensions of the camera image
            w = obj.snapDim(1);
            h = obj.snapDim(2);
            
            % Particles pixel position and yaw angle
            xStartCell = num2cell(round(partImgCenterpxPos(:,1) - w/2));
            yStartCell = num2cell(round(partImgCenterpxPos(:,2) - h/2));

            W = obj.AIM.mapDim(1);
            H = obj.AIM.mapDim(2);
            
            % Cropping and rotating each of particles with using 'cellfun'
            % rather than a 'for loop'. Resulting array is a 'cell array'.
            PartImage = cellfun(@(x,y) imrotate(obj.AIM.I(max(min(y:y+h-1,H),1), max(min(x:x+w-1,W),1), :), yaw, 'crop'),xStartCell,yStartCell,'UniformOutput',false);

        end

        function inlierIdx = ImageMatching(obj,UAVimage,PartImage)
            %%This function return matched feature points between UAV Image
            %%and Particle Image after remove outliers. Different feature
            %%extraction and feature matching algorithm can be used. For
            %%now, this function uses 'ORB' feature extractor and 'kNN'
            %%based feature matching method. In the next version learning
            %%based feature extractor 'Superpoint' and feature matcher
            %%'LightGlue' will be use. 'inlierIdx' is a cell array
            %%contais boolean(0 or 1) array in each elements that
            %%demonstrate mathced feature inlier or not. We then count how
            %%many inliers points particles have for finding likelihood.

            %Feature extraction of UAV image
            [featuresUAV,validpointsUAV] = extractFeatures(UAVimage,detectORBFeatures(UAVimage,"NumLevels",obj.num_level));

            % tic
            % Feature extraction of particles images using 'cellfun' rather
            % then 'forloop'
            [featuresPart,validpointsPart] = cellfun(@(x) extractFeatures(x,detectORBFeatures(x,"NumLevels",obj.num_level)),PartImage,'UniformOutput',false);
            % toc

            % tic
            % Featuer matching between UAV and Particles images
            indexPairs = cellfun(@(x) matchFeatures(featuresUAV,x),featuresPart,'UniformOutput',false,'ErrorHandler',@obj.errorFunc);
            % toc

            % tic
            % Outlier removal prosedure
            [~,inlierIdx,~] = cellfun(@(PartValid,idxPair) estgeotform2d(PartValid(idxPair(:,2),:),validpointsUAV(idxPair(:,1),:),"similarity"),validpointsPart,indexPairs,'UniformOutput',false,'ErrorHandler',@obj.errorFunc);
            % toc
        end

        function [B, A, C] = errorFunc(varargin)
            %This is a helper function for preventing 'cellfun' error
            % massage occurs
            A = 0; 
            B = 0;
            C = 0;
        end
       
    end
end
