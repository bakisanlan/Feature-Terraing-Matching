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

    end

    methods

        function obj = DatabaseScanner(varargin)
            %DatabaseScanner Construct an instance of this class

            % Call parent class constructor
            obj.slidingTypeInterp = true;

        end

        function [Zs, Zw]= find_likelihood(obj,orientation,particles,UAVImage)
            % FETCHPARTMEASUREMENT function will be called by State Estimators(PSO,PF)
            % to get measurement of particles. Function outputs are 2D
            % array that's size are NxM. N is number of particles and M is
            % number of point in point cloud measurement.

            % Assume altitude and orientation are known.
            obj.positionLiDAR(3) = alt;
            obj.orientationLiDAR = orientation;

            Zr_s = ptCloudRadar.s.Location(:,3);
            Zr_w = ptCloudRadar.w.Location(:,3);

            numParticles = numel(particles(:,1));
            numPointPC   = numel(Zr_w);

            Zs = zeros(numParticles,numPointPC);
            Zw = zeros(numParticles,numPointPC);

            % Raycast or Sliding methods will be used depending on raycastFlag
            if raycastFlag
                for i=1:numParticles

                    xParticle = particles(i,1);
                    yParticle = particles(i,2);
                    obj.positionLiDAR([1, 2]) = [xParticle yParticle];

                    % Point cloud(LIDAR like) or single point(Altimeter
                    % like) measurement by each particle.
                    if length(Zr_w) == 1
                        obj.scanTerrainOneRay(false)
                    else
                        obj.scanTerrain(false)
                    end

                    Zs(i,:) = obj.ptCloud.s.Location(:,3);
                    Zw(i,:) = obj.ptCloud.w.Location(:,3);
                end

            else
                PCSensor = ptCloudRadar.s.Location;
                if obj.slidingTypeInterp
                    [Zs, Zw] = obj.slidingInterp(particles,PCSensor);
                else
                    [Zs, Zw]= obj.slidingOrthogRay(particles,PCSensor);
                end
            end

        end

        function snappedImage = snapImage(obj,imgCenterWorldPos,ROT)

            imgCenterpxPos = imgCenterWorldPos * (1/obj.AIM.mp);


            


            
        end

        function scanTerrainOneRay(obj,flagPlot)
            %SCANTERRAIN shoot single ray to terrain in the
            %direction of -z axis vector of sensor frame. This emulates
            %banked Radar Altimeter. This function is identical with
            %terrain.RadarAltimeter3DMesh.scanTerrain

            [xs,ys,zs,xw,yw,zw]      = raycast(obj, 90, 0);

            obj.ptCloud.s = pointCloud([xs ys zs]);
            obj.ptCloud.w = pointCloud([xw yw zw]);

            if flagPlot && ~isempty(obj.hFigure)
                figure(obj.hFigure);
                plot3(Xw,Yw,Zw,'r.','MarkerSize',10,'DisplayName','LiDAR Scans');
            end
        end

        function [zs, zw] = slidingOrthogRay(obj,particles,PCSensor)
            %SLIDINGORTHOGRAY takes inputs as particles and point cloud
            %measurement in sensor frame then finds zw 2D array with single
            %point orthogonal ray.

            %Function outputs are 2D array that's size are NxM. N is number
            % of particles and M is number of point in point cloud measurement.

            % Save the sensor location and range
            originalPos     = obj.positionLiDAR;
            originalOrient  = obj.orientationLiDAR;
            originalRange   = obj.rayRange;

            [num_Part, ~] = size(particles);
            numPointPCsensor = length(PCSensor(:,1));

            zw = zeros(num_Part,numPointPCsensor);
            zs = zeros(num_Part,numPointPCsensor);

            psi_s = obj.orientationLiDAR(1);  % Yaw angle
            phi_s = obj.orientationLiDAR(3);  % Banking angle

            % Iterating for every particles and for every point in
            % point cloud.
            for i=1:num_Part
                for j=1:numPointPCsensor

                    obj.positionLiDAR = particles(i,:)';
                    Point_xyz_s = PCSensor(j,:)'; %jth point in sensor frame

                    % Converting jth point location to world frame
                    Point_xyz_w = terrain.AbstractRayCasting3D.wTs(obj.positionLiDAR) * ...
                        terrain.AbstractRayCasting3D.wRs(psi_s,phi_s,true) * [Point_xyz_s;1];

                    % Move to radar point on z direction of world ENU to aircraft
                    % level
                    Point_xyz_w(3) = obj.positionLiDAR(3);

                    % Move the sensor to the desired point and ray cast from there
                    obj.positionLiDAR   = Point_xyz_w(1:3);
                    obj.orientationLiDAR(3) = 0;   % Ray cast orthogonally without a bank angle
                    obj.rayRange        =  obj.positionLiDAR(3);  % limits terrain area to a minimum
                    [~,~,zs(i,j),~,~,zw(i,j)]      = raycast(obj, 90, 0);

                    % Restore the sensor orientation and range
                    obj.orientationLiDAR = originalOrient;
                    obj.rayRange         = originalRange;
                end
            end
            % Restore the sensor location
            obj.positionLiDAR   = originalPos;
        end

        function [zs, zw] = slidingInterp(obj,particles,PCSensor)
            % SLIDINGINTERP takes inputs as particles and point cloud
            %measurement in sensor frame then finds zw 2D array without
            %raycasting through interpolation method. This method do not
            %use any for loop as like slidingOrthogRay.

            %Function outputs are 2D array that's size are NxM. N is number
            % of particles and M is number of point in point cloud measurement.

            psi_s = obj.orientationLiDAR(1);
            phi_s = obj.orientationLiDAR(3);


            [num_Part, col] = size(particles);
            numPointPCsensor = length(PCSensor(:,1));

            %Below codes, provide to get particle measurement without 
            %using any 'for loop' through cell and matrix operations.

            %Nx1 cell array that consatins states of each particles
            particlesPosWorldCell = mat2cell(particles,ones(1,num_Part), col);
            
            % Nx1 cell array that contains augmented translation matrix(4x4)
            transMatCell = cellfun(@(x) terrain.AbstractRayCasting3D.wTs(x),particlesPosWorldCell,'UniformOutput',false);
            
            % Augmented rotation matrix(4x4)
            [~, rotMat] = terrain.AbstractRayCasting3D.wPOSEs([0;0;0],[0 0 0],psi_s,phi_s);
            
            % Nx1 cell array that contains same augmented rotation matrix(4x4) in each cell.
            rotMatCell = repmat({rotMat}, num_Part, 1);
            
            %Force PCSensor to be augmented form Mx4
            PCSensor = [PCSensor ones(numPointPCsensor,1)]; 
            
            %To find point cloud position of particles in world frame,
            %rotation and translation is made by 'cellfun'.
            particlesPCWorldCell = cellfun(@(x,y) ((x*y)*PCSensor')',transMatCell,rotMatCell,'UniformOutput',false);
            
            %Convert cell array to 3D array
            particlesPCWorldArray = cat(3,particlesPCWorldCell{:});
            
            %Remove augmented columns from particlesPCWorldArray
            particlesPCWorldArray(:,4,:) = [];
            
            %2D gridded interpolation with DTED grid map.
            zw = interp2(obj.DTED{1},obj.DTED{2},obj.DTED{3},particlesPCWorldArray(:,1,:),particlesPCWorldArray(:,2,:)); %Mx3xN
            particlesPCWorldArray(:,3,:) = zw;
            %Converting zw 3D(Mx1xN) array to 2D(NxM)
            zw = reshape(zw,[numPointPCsensor, num_Part])';

            % RESOLVE LATER: Zs is necessary?
            zs = NaN;
        end
    end

end
