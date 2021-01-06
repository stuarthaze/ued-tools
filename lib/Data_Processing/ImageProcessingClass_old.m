classdef ImageProcessingClass
    
    properties
        IMG_raw
        IMG
        useIntensityCalibration
        useCropping
        
        useMask
        useRadialBG
        processingOptions
        IMG_processed
        IMG_cropped
        IMG_bgsub
        IMG_filtered
        BG
        backgroundOptions
        filter
        cropping
        calibration
        binning
        beamstop
        beamstopOptions
        mask
        maskOptions
        imgsizeFullFrame
        imgsize
        centreFullFrame
        centre
        Rmatrix
        peaks
        peaksAuto
        displayOptions
    end %properties
    
    methods
        %Constructor
        function obj = ImageProcessingClass(IMGin)
            obj.IMG_raw = IMGin;
        end
        %------------
        % Read data
        function obj = readIMG(obj,dirData,string2match)
            obj.IMG_rawData = readSPE_images(dirData,string2match);
            if ndims(obj.IMG_rawData) == 3
                numImages = size(obj.IMG_rawData,3);
                disp(['Warning: ',num2str(numImages),'images found']);
            end
        end
        function obj = readIMG_fromBackups(obj,backupDir,string2match)
            threshold = 3;
            ratioNormalQuantile = 0.5;
            displayData = 0;
            displayInfo = 0;
            [obj.IMG_rawData, obj.readStats.numPixelsContributingIMG] = processImagesInDirectory_16bit_Hybrid(...
                backupDir,string2match,threshold,ratioNormalQuantile,displayData,displayInfo);
        end
        function obj = readBG(obj,dirBG,string2match)
            obj.IMG_detctorBG = readSPE_images(dirBG,string2match);
            if ndims(obj.IMG_detctorBG) == 3
                disp(['Warning: ',num2str(numImages),'images found']);
            end
        end
        function obj = readBG_fromBackups(obj,backupDir,string2match)
            threshold = 3;
            ratioNormalQuantile = 0.5;
            displayData = 0;
            displayInfo = 0;
            [obj.IMG_detectorBG, obj.readStats.numPixelsContributingBG] = processImagesInDirectory_16bit_Hybrid(...
                backupDir,string2match,threshold,ratioNormalQuantile,displayData,displayInfo);
        end
            
        %-------------
        % Image processing
        function obj = processRawData(obj)
            obj.imgsizeFullFrame = size(obj.IMG_rawData);
            IMG = obj.IMG_rawData;
            if obj.processingOptions.subtractDetectorBG
                IMG = IMG - obj.IMG_detectorBG;
            end
            if obj.processingOptions.useIntensityCalibration
                obj.calibration.intensities.fullIMG = createDetectorCalibrationMat(obj.imgsizeFullFrame,obj.binning,obj.processingOptions.nameCalibrationMat);
                IMG = IMG.*obj.calibration.intensities.fullIMG;
            end
            obj.IMG_processed = IMG;
        end
        %------------------------------------
        % Cropping
        function obj = setCroppingDetials(obj)
            obj.cropping.x = obj.centreFullFrame.x + (1:obj.imgsize.x)- round(obj.imgsize.x/2);
            obj.cropping.y = obj.centreFullFrame.y + (1:obj.imgsize.y)- round(obj.imgsize.y/2);
            obj.centre.x = obj.centreFullFrame.x - obj.cropping.x(1) + 1;
            obj.centre.y = obj.centreFullFrame.y - obj.cropping.y(1) + 1;
        end
        function obj = cropIMG(obj)
            obj.IMG_cropped = obj.IMG_processed(obj.cropping.y, obj.cropping.x);
        end
        %------------------------------------------
        % Mask
        function obj = createBeamStopMask(obj)
            obj.beamstop = findBeamBlock(obj.IMG_cropped,obj.centre,obj.beamstopOptions.nPixelExpand);
        end
        
        function obj = createMask(obj)
            obj.mask = zeros([obj.imgsize.y,obj.imgsize.x]);
            if obj.maskOptions.useBeamstop && ~isempty(obj.beamstop)
                obj.mask = obj.mask + obj.beamstop;
            end
        end
        %-------------------------------
        % Radial BG
        function obj = calcRadialBG(obj)
            [obj.BG.radial2d, obj.BG.radial1d] = radialQuantile_BG(obj.IMG_cropped,...
                                                obj.centre,obj.mask,obj.options.quantileRadialBG);
        end
        %-----------------------------
        function obj = calcRmatrix(obj)
            X = 1:obj.imgsize.x;
            Y = 1:obj.imgSize.y;
            [XX,YY] = meshgrid(X,Y);
            RR2rd = (XX-obj.centre.x).^2 + (YY-obj.centre.y).^2;
            obj.Rmatrix = sqrt(RR2rd);
        end
        %-------------
        % Display data
        function display_IMG_processed(obj,options)
            obj.displayImage(obj.IMG_processed,options);
        end
        function display_IMG_cropped(obj,options)
            obj.displayImage(obj.IMG_cropped,options);
        end

        function displayImage(obj,img,options)
            if options.useMask
                img = img.*(1-obj.mask);
            end
            if options.scalePower
                img(img < 0) = 0;
                img = img.^options.scalePower;
            end
            scale = options.scale;
            if scale == 0
                scale = max(max(img));
            end
            if length(scale) == 1
                scalemin = 0;
                scalemax = scale;
            elseif length(scale) == 2
                scalemin = scale(1);
                scalemax = scale(2);
            end

            figure();
            if options.invert
                imagesc(-img,[-scalemax,-scalemin]); 
            else
                imagesc(img,[scalemin,scalemax]); 
            end
            pbaspect([size(img,2),size(img,1),1]);
            colormap(options.colormap);
        end
    end %methods
end %class
%----------------
% Utility Functions

