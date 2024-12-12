classdef AerialImageModel < handle
    %DigitalElevationModel Manages a Digital Elevation Model (DEM) for
    %terrain scan and matching tasks.
    %   Processes a DTED file gathered by NASA SRTM.

    properties
        % Map dimensions
        mapDim;

        % A and R
        I  % base image in mxnxc sized
        mp % meter/pixel ratio
    end

    methods
        function obj = AerialImageModel(area)
            %DigitalElevationModel Constructs a DEM terrain object.
            
            % Load elevation data from fill
            if strcmpi(area,'ITU')
                % DTED at (41N, 29E) of the Bosphorus.
                obj.loadData(fullfile(fileparts(mfilename('fullpath')),...
                    'data/itu_map_square.jpg'));
            elseif strcmpi(area,'GC')
               obj.loadData(fullfile(fileparts(mfilename('fullpath')),...
                    'data/n36_w113_1arc_v3.dt2'));
            elseif strcmpi(area, 'MR')
                % load Marmara Region with 40-42 latitudes and 27-30 longitudes
                lats = 40:41; % DTED latitudes range
                lons = 27:29; % DTED longitudes range
                obj.loadBatchData(lats, lons);
            else
                error('Enter a valid area name.')
            end
        end

        function loadData(obj,filename)
            %loadData Load data from file

            obj.I  = imread(filename);
            obj.mapDim = size(obj.I);
            obj.mp = 725/obj.mapDim(1);

        end

        function loadBatchData(obj, lats, lons)
            % load batch of latitudes longitudes
           
        end

        function Icropped = slice(obj,xmin,ymin,width,height)
            %slice Return a slice of the loaded DTED data.
            Icropped = imcrop(obj.I,[xmin,ymin,width,height]);
            
        end


        function visualizeBaseAerialImage(obj)
            %visualizeDTED visualizes loaded DTED file.

            figure; clf;
            imshow(obj.I)

        end

        function visualizeDTEDSlice(obj,xmin,ymin,width,height)
            %visualizeDTED visualizes a given slice of the loaded DTED
            %file.

            figure; clf;
            Icropped = slice(obj,xmin,ymin,width,height);
            imshow(Icropped)
        end

        function visualizeDTEDSliceGray(obj,xmin,ymin,width,height)
            %visualizeDTED visualizes a given slice of the loaded DTED
            %file in Gray color.

            figure; clf;
            Icropped = slice(obj,xmin,ymin,width,height);
            Igray    = im2gray(Icropped);
            imshow(Igray)
        end

    end
end