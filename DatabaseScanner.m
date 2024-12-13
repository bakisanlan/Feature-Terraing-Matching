classdef DatabaseScanner < handle
    %DATABASESCANNER Class wil be used by State Estimators(PSO,PF) for
    %getting measurements of particles through scannig the offline
    %database map(DTED).

    % Radar pattern of aircraft can be point cloud measurement data by LIDAR
    % like sensor or single point measurement data by RadarAltimeter sensor
    % like.

    % Measurement of particles can be found by below two methods those are
    % raycasting, sliding.
    %       - Each particle raycasts from its current reference map
    %         position (xw,yw) to get point cloud (Xp,Yp,Zp) represented in
    %         the particle's frame. Elevations Zp are then compared with
    %         Zm. Currently, the (Xp,Yp) information is not compared with
    %         (Xm,Ym).
    %       - Each particle uses the (Xm,Ym) pattern and its current
    %         reference map position (xw,yw) to directly read Zp elevation
    %         values at (xw,yw)+(Xm,Ym) to get a point cloud (Xm,Ym,Zp).
    %         The Zp elevations are then compared with Zm. Thus no
    %         raycasting is performed.

   

    properties
        % There are two different type of sliding methods those are 
        % interpolation method, single orthogonal raycasting. TODO: second
        % one can be deleted becasue of interp method is very fast.
        snapDim
        AIM
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
            % FETCHPARTMEASUREMENT function will be called by State Estimators(PSO,PF)
            % to get measurement of particles. Function outputs are 2D
            % array that's size are NxM. N is number of particles and M is
            % number of point in point cloud measurement.

            PartImage = obj.snapPartImage(partImgCenterWorldPos,yaw);

            inlierIdx = obj.ImageMatching(UAVImage,PartImage);

            numMatchedFeature = cellfun(@(x) sum(x),inlierIdx);

            % deal later of likelihood function
            partlikelihood = numMatchedFeature;
        end

        function UAVImage = snapUAVImage(obj,UAVWorldPos,yaw)
        
            % rescale y component due to direction of y axis of image frame
            UAVWorldPos(2) = -UAVWorldPos(2) + 200;

            UAVImgCenterpxPos = round(UAVWorldPos * (1/obj.AIM.mp));

            % Dimensions of the camera image
            w = obj.snapDim(1);
            h = obj.snapDim(2);

            xStart = UAVImgCenterpxPos(1) - w/2;
            yStart = UAVImgCenterpxPos(2) - h/2;

            UAVImage = imrotate(obj.AIM.I(yStart:yStart+h-1, xStart:xStart+w-1, :), yaw, 'crop');
    
        end

        function PartImage = snapPartImage(obj,partWorldPos,yaw)

            % rescale y component due to direction of y axis of image frame
            partWorldPos(:,2) = -partWorldPos(:,2) + 200;

            % Convert (N,E)m dimension to (u,v)px
            partImgCenterpxPos = round(partWorldPos * (1/obj.AIM.mp));

            % Dimensions of the camera image
            w = obj.snapDim(1);
            h = obj.snapDim(2);
            
            % UAV position and yaw angle
            xStartCell = num2cell(round(partImgCenterpxPos(:,1) - w/2));
            yStartCell = num2cell(round(partImgCenterpxPos(:,2) - h/2));
            
            % tic
            PartImage = cellfun(@(x,y) imrotate(obj.AIM.I(y:y+h-1, x:x+w-1, :), yaw, 'crop'),xStartCell,yStartCell,'UniformOutput',false);
            % toc

        end

        function inlierIdx = ImageMatching(obj,UAVimage,PartImage)

            %feature extraction
            
            [featuresUAV,validpointsUAV] = extractFeatures(UAVimage,detectORBFeatures(UAVimage,"NumLevels",obj.num_level));

            % tic
            [featuresPart,validpointsPart] = cellfun(@(x) extractFeatures(x,detectORBFeatures(x,"NumLevels",obj.num_level)),PartImage,'UniformOutput',false);
            % toc

            % tic
            indexPairs = cellfun(@(x) matchFeatures(featuresUAV,x),featuresPart,'UniformOutput',false);
            % toc

            % tic
            [~,inlierIdx,~] = cellfun(@(PartValid,idxPair) estgeotform2d(PartValid(idxPair(:,2),:),validpointsUAV(idxPair(:,1),:),"similarity"),validpointsPart,indexPairs,'UniformOutput',false);
            % toc
        end

       
    end
end
