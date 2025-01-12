classdef AerialImageModel < handle
    %DigitalElevationModel Manages a Digital Elevation Model (DEM) for
    %terrain scan and matching tasks.
    %   Processes a DTED file gathered by NASA SRTM.

    properties

        mapDim;         % Map dimensions
        num_level       % Orb feature descriptor parameter
        I               % base image in mxnxc sized
        mp              % meter/pixel ratio
        featuresBase    %Feature descriptor of Base Image 
        validpointsBase %Interest point of Base Image
    end

    methods
        function obj = AerialImageModel(area)            
            % Load elevation data from fill
            if strcmpi(area,'ITU')
                % ITU Location Satellite Image
                obj.loadData(fullfile(fileparts(mfilename('fullpath')),...
                    'data/itu_map_square.jpg'));
            else
                error('Enter a valid area name.')
            end
        end

        function loadData(obj,filename)
            %loadData Load data from file

            obj.I  = im2double(im2gray(imread(filename))); %gray scale double type image array
            obj.mapDim = size(obj.I);
            obj.mp = 725/obj.mapDim(1);     % 725meters ~= mapDim(1), meter/pixel ratio
            obj.num_level = 1;
            [obj.featuresBase,obj.validpointsBase] = extractFeatures(obj.I,detectORBFeatures(obj.I,"NumLevels",obj.num_level));


        end

        function Icropped = slice(obj,xmin,ymin,width,height)
            %slice Return a slice of the loaded IMAGE data.
            Icropped = imcrop(obj.I,[xmin,ymin,width,height]);
            
        end


        function visualizeBaseAerialImage(obj)
            %visualizeDTED visualizes loaded IMAGE file.

            figure; clf;
            imshow(obj.I)

        end

        function visualizeSlicedAerialImage(obj,xmin,ymin,width,height)
            %visualizeSlicedAerialImage visualizes a given slice of the
            %loaded IMAGE file.

            figure; clf;
            Icropped = slice(obj,xmin,ymin,width,height);
            imshow(Icropped)
        end

        function visualizeSlicedGrayAerialImage(obj,xmin,ymin,width,height)
            %visualizeSlicedGrayAerialImage visualizes a given slice of the
            %loaded IMAGE file in Gray color.

            figure; clf;
            Icropped = slice(obj,xmin,ymin,width,height);
            Igray    = im2gray(Icropped);
            imshow(Igray)
        end

    end
end