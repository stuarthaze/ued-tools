classdef ImageProcessingClass2
    
    properties
        IMG
        IMG_raw
        IMG_unfiltered
        % Processing Options
        useIntensityCalibration
        useCropping
        useZeroThreshold
        maskBeamstop4bg
        maskRmin
        maskRmax
        useRadialBG
        useFilter
        
        BG_radial2d
        BG_radial1d
        quantileRadialBG
        filter
        cropping
        cropsize
        calibrationName
        calibrationIMG
        binning
        beamstop
        beamstopEdges
        beamstopExpansion
        mask_bgSubtraction
        mask_peakFinding
        imgsizeRaw
        imgsize
        centreRaw
        centre
        Rmatrix
        peaks
        numPeaks
        peakDiam
        peaksAuto
        ROIstack
        ROIvalues
        ROI_numPixels
        
    end %properties
    
    methods
        %Constructor
        function obj = ImageProcessingClass2(IMGin)
            obj.IMG_raw = IMGin;
            obj.IMG = obj.IMG_raw;
            obj.imgsizeRaw = size(obj.IMG_raw);
            obj.useIntensityCalibration = 0;
            obj.useCropping = 0;
            obj.maskBeamstop4bg = 0;
            obj.useRadialBG = 0;
        end
        %-------------------------------
        % Processing Functions
        function obj = applyDetectorCalibration(obj)
            if isempty(obj.calibrationIMG)
                obj = obj.loadDetectorCalibration();
            end
            obj.IMG = obj.IMG_raw.*obj.calibrationIMG;
        end
        
        % Process Data Automatically
        function obj = processIMG(obj)
            obj.IMG = obj.IMG_raw;
            %-- Calibration
            if obj.useIntensityCalibration
                obj = obj.applyDetectorCalibration();
            end
            %-- Cropping
            if obj.useCropping
                obj = obj.cropImage();
            else
                obj.centre = obj.centreRaw;
            end
%             obj.imgsize = size(obj.IMG);
            
            %-- Radial BG
            if obj.useRadialBG
                obj = obj.calcRadialBG();
                obj = obj.subtractRadialBG();
            end
            if obj.useZeroThreshold
                obj.IMG(obj.IMG < 0) = 0;
            end
        end
        %----------------------------------
        function obj = applyZeroThreshold(obj)
            obj.IMG(obj.IMG < 0) = 0;
        end
        %----------------------------------
        % Detector Calibration
        function obj = loadDetectorCalibration(obj)
            calibrationDet = createDetectorCalibrationMat(obj.imgsizeRaw,obj.binning,obj.calibrationName);
            calibrationMedian = median(median(calibrationDet,'omitnan'));
            obj.calibrationIMG = medfilt2(calibrationDet/calibrationMedian);
        end
        %----------------------------------
        % Peak Finding
        function obj = findPeaks(obj,options)
            PKS = findPeaks2d(obj.IMG,options.windowsize);
            mask = zeros(size(obj.IMG));
            if options.maskBeamstop
                if isempty(obj.beamstop)
                    obj = obj.createBeamstopMask();
                end                    
                Kern = createCircle(options.minDistToBeamStop);
                mask = conv2(double(obj.beamstop),Kern,'same');
            end
            if options.Rmin
                mask(obj.Rmatrix < options.Rmin) = 1;
            end
            if options.Rmax
                mask(obj.Rmatrix > options.Rmax) = 1;
            end
            mask(mask > 1) = 1;
            obj.mask_peakFinding = mask;
            PKS_masked = PKS.*(1-mask);
            [peakHeights, peakPositionsVec] = sortPeaks(PKS_masked);
            nPeaksFound = sum(peakHeights > 0);
            if options.numPeaksAuto < nPeaksFound
                nPeaks = options.numPeaksAuto;
            else
                nPeaks = nPeaksFound;
            end
            obj.peaksAuto.x = ceil(peakPositionsVec(1:nPeaks)/size(obj.IMG,1));
            obj.peaksAuto.y = mod(peakPositionsVec(1:nPeaks),size(obj.IMG,1));
%             obj.peaksAuto.vals = peakHeights(1:nPeaks);
            obj.peaks = obj.peaksAuto;
            obj.numPeaks = nPeaks;
            obj.peakDiam = options.diameter;
        end
        
        function obj = readPeaksFromAuto(obj)
            obj.peaks = obj.peaksAuto;
            obj.numPeaks = length(obj.peaks.vals);
        end
        function obj = deletePeaks(obj,peakIndices)
            peaks2keep = 1:obj.numPeaks > 0;
            peaks2keep(peakIndices) = 0;
            obj.peaks.x = obj.peaks.x(peaks2keep);
            obj.peaks.y = obj.peaks.y(peaks2keep);
%             obj.peaks.vals = obj.peaks.vals(peaks2keep);
            obj.numPeaks = length(obj.peaks.x);
        end
        function obj = addPeaks(obj,Xnew,Ynew)
            obj.peaks.x = [obj.peaks.x, Xnew];
            obj.peaks.y = [obj.peaks.y, Ynew];
            obj.numPeaks = length(obj.peaks.x);
        end
        %------------------------------------
        % Cropping
        function obj = cropImage(obj)
            if isempty(obj.cropsize) && isempty(obj.cropping)
                disp('Error: Cropping details missing');
                return
            elseif isempty(obj.cropsize)
                obj.cropsize(2) = length(obj.cropping.x);
                obj.cropsize(1) = length(obj.cropping.y);
                obj.imgsize.x = obj.cropsize(2);
                obj.imgsize.y = obj.cropsize(1);
            elseif isempty(obj.cropping)
                obj = obj.updateCroppingVals();
            end
                
            obj.centre.x = obj.centreRaw.x - obj.cropping.x(1) + 1;
            obj.centre.y = obj.centreRaw.y - obj.cropping.y(1) + 1;
            obj.IMG = obj.IMG(obj.cropping.y,obj.cropping.x);
            obj = obj.calcRmatrix();
        end
        
        function obj = updateCroppingVals(obj)
            obj.imgsize.x = obj.cropsize(2);
            obj.imgsize.y = obj.cropsize(1);
            obj.cropping.x = obj.centreRaw.x + (1:obj.imgsize.x)- round(obj.imgsize.x/2);
            obj.cropping.y = obj.centreRaw.y + (1:obj.imgsize.y)- round(obj.imgsize.y/2);
        end
        %--------------------------------------------
        
        function obj = setCentre(x,y)
            obj.centre.x = x;
            obj.centre.y = y;
            obj = obj.calcRmatrix();
        end
        %------------------------------------------
        % Mask
        function obj = createBeamstopMask(obj)
            [obj.beamstop, obj.beamstopEdges] = findBeamBlock(obj.IMG,obj.centre,obj.beamstopExpansion);
        end
%         
%         function obj = createMask(obj)
%             obj.mask = zeros([obj.imgsize.y,obj.imgsize.x]);
%             if obj.maskOptions.useBeamstop && ~isempty(obj.beamstop)
%                 obj.mask = obj.mask + obj.beamstop;
%             end
%         end
        %-------------------------------
        % Radial BG
        function obj = calcRadialBG(obj)
            obj.mask_bgSubtraction = zeros(size(obj.IMG)) > 1; %Creates a logical array of 1's
            if obj.maskBeamstop4bg
                if isempty(obj.beamstop)
                    obj = obj.createBeamstopMask();
                end
                obj.mask_bgSubtraction = obj.beamstop;
            end
            obj = obj.calcRmatrix();
            if obj.maskRmin
                obj.mask_bgSubtraction = obj.mask_bgSubtraction | (obj.Rmatrix <= obj.maskRmin);
            end
            if obj.maskRmax
                obj.mask_bgSubtraction = obj.mask_bgSubtraction | (obj.Rmatrix >= obj.maskRmax);
            end

            [obj.BG_radial2d, obj.BG_radial1d] = radialQuantile_BG(obj.IMG,...
                                                obj.centre,obj.mask_bgSubtraction,obj.quantileRadialBG);
        end
        
        function obj = subtractRadialBG(obj)
            obj = obj.calcRadialBG();
            obj.IMG = obj.IMG - obj.BG_radial2d;
        end
        %-----------------------------
        function obj = calcRmatrix(obj)
            X = 1:size(obj.IMG,2);
            Y = 1:size(obj.IMG,1);
            [XX,YY] = meshgrid(X,Y);
            RR2rd = (XX-obj.centre.x).^2 + (YY-obj.centre.y).^2;
            obj.Rmatrix = sqrt(RR2rd);
        end
        %---------------------------
        % Filtering
        function obj = filterIMG(obj,filterOpt)
            if ~obj.useFilter
                if isempty(obj.IMG_unfiltered)
                    return
                else
                    obj.IMG = obj.IMG_unfiltered;
                    return
                end
            end
            if isempty(obj.IMG_unfiltered)
                obj.IMG_unfiltered = obj.IMG;
            elseif filterOpt.resetFilter
                obj.IMG = obj.IMG_unfiltered;
            end
            if filterOpt.median
                if length(filterOpt.median)==1
                    MxN = [filterOpt.median, filterOpt.median];
                else
                    MxN = filterOpt.median;
                end
                obj.IMG = medfilt2(obj.IMG,MxN);
            end
            if filterOpt.nearestNeighbour
                obj.IMG = fNearestNeighbourFilter(obj.IMG,filterOpt.nearestNeighbour);
            end
            if filterOpt.outlierThreshold
                obj.IMG = outlierFilter2dGauss(obj.IMG,filterOpt.outlierThreshold);
            end

            if filterOpt.sigmaGaussian
                obj.IMG = imgaussfilt(obj.IMG,filterOpt.sigmaGaussian);
            end
            if filterOpt.sigmaMexhat
                IMGmex = fMexicanHatFilter2d(obj.IMG,filterOpt.sigmaMexhat);
                IMGmex = IMGmex/max(max(IMGmex));
                IMGmex(IMGmex < 0) = 0;
                IMGnorm = obj.IMG/max(max(obj.IMG));
                obj.IMG = filterOpt.ratioMexhat*IMGmex + (1-filterOpt.ratioMexhat)*IMGnorm;
            end
            if filterOpt.sigmaGaussian2
                obj.IMG = imgaussfilt(obj.IMG,filterOpt.sigmaGaussian2);
            end
            if filterOpt.radialExponent
                obj.IMG = obj.IMG.*(obj.Rmatrix.^filterOpt.radialExponent);
            end
        end
        %-------------
        % Display data
        function displayIMG(obj,displayOptions)
            img = obj.IMG;
            if displayOptions.maskBeamstop
                img = img.*(1-double(obj.beamstop));
            end
            if displayOptions.gamma
                img(img < 0) = 0;
                img = img.^displayOptions.gamma;
            end
            scale = displayOptions.scale;
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
            imagesc(img,[scalemin,scalemax]); 
            cmap = colormap(displayOptions.colormap);
            if displayOptions.invert
                colormap(flipud(cmap));
            end

            pbaspect([size(img,2),size(img,1),1]);
            

%             figure();
%             if displayOptions.invert
%                 imagesc(-img,[-scalemax,-scalemin]); 
%             else
%                 imagesc(img,[scalemin,scalemax]); 
%             end
%             pbaspect([size(img,2),size(img,1),1]);
%             colormap(displayOptions.colormap);
            if displayOptions.viewPeaks
                peakRadius = obj.peakDiam*0.5;
                if ~isempty(displayOptions.peaksToView)
                    peakNumbers = displayOptions.peaksToView;
                else
                    peakNumbers = 1:obj.numPeaks;
                end
                xpos = obj.peaks.x(peakNumbers);
                ypos = obj.peaks.y(peakNumbers);
                addCirclesToPlot(xpos, ypos, peakRadius, displayOptions.peakTextSize, displayOptions.peakTextColor, peakNumbers);
            end
        end
        %----------------------
        % Generate ROI masks for data extraction
        function obj = generateROImasks(obj,diam)
            obj.ROIstack = zeros([obj.imgsize.y,obj.imgsize.x,obj.numPeaks]);
            for ipeak = 1:obj.numPeaks
                ROI = zeros(size(obj.IMG));
                ROI(obj.peaks.y(ipeak),obj.peaks.x(ipeak)) = 1;
                obj.ROIstack(:,:,ipeak) = convolveCircle(ROI,diam);
            end
            obj.ROI_numPixels = sum(sum(obj.ROIstack(:,:,1)));
        end
        % Calculate ROI values
        function obj = calcROIvalues(obj)
            obj.ROIvalues = zeros([1,obj.numPeaks]);
            if isempty(obj.IMG_unfiltered)
                img = obj.IMG;
            else
                img = obj.IMG_unfiltered;
            end
            for ii = 1:obj.numPeaks
                obj.ROIvalues(ii) = sum(sum(img.*obj.ROIstack(:,:,ii)))/obj.ROI_numPixels;
            end
        end
        %--------------------------------
        % Save parameters to text files
        function savePeaks(obj)
            newFolder = 'Peak Finding';
            if ~isdir(['.\',newFolder])
                mkdir(newFolder);
            end
            if isempty(obj.IMG_unfiltered)
                img = obj.IMG;
            else
                img = obj.IMG_unfiltered;
            end
            X_Y_Intensity = [obj.peaks.x',obj.peaks.y',obj.ROIvalues'];
            save(['./',newFolder,'/X_Y_Intensity.txt'],'X_Y_Intensity','-ascii');
            XY_centre = [obj.centre.x, obj.centre.y];
            save(['./',newFolder,'/XY_centre.txt'],'XY_centre','-ascii');
            X_cropping_min_max = [obj.cropping.x(1), obj.cropping.x(end)];
            Y_cropping_min_max = [obj.cropping.y(1), obj.cropping.y(end)];
            save(['./',newFolder,'/X_cropping_min_max.txt'],'X_cropping_min_max','-ascii');
            save(['./',newFolder,'/Y_cropping_min_max.txt'],'Y_cropping_min_max','-ascii');
            save(['./',newFolder,'/IMG_cropped.mat'],'img');
        end
    end %methods
end %class
%----------------
% Utility Functions
function C = createCircle(radius)
    X = -radius:radius;
    [XX,YY] = meshgrid(X,X);
    RR = sqrt(XX.^2 + YY.^2);
    C = RR <= radius;
    C = double(C);
end

function addCirclesToPlot(xCen, yCen, radii, textSize, textColor, numbering)
    if nargin > 3
        labelcircles = 1;
    else
        labelcircles = 0;
    end
    nx = length(xCen);
    ny = length(yCen);
    nr = length(radii);
    if nx ~= ny
        disp('error in centres');
        return
    end
    if nx == nr
        xCenVec = xCen;
        yCenVec = yCen;
        rVec = radii;
    elseif nr == 1
        xCenVec = xCen;
        yCenVec = yCen;
        rVec = radii*ones([1,nx]);
    elseif nx == 1
        xCenVec = xCen*ones([1,nr]);
        yCenVec = yCen*ones([1,nr]);
        rVec = radii;
    end
    nCircles = length(xCenVec);
    th = (0:0.01:1)*2*pi;
    x_circle = cos(th);
    y_circle = sin(th);
    %-------
    hold on;
    for ii = 1:nCircles
        x = xCenVec(ii) + rVec(ii)*x_circle;
        y = yCenVec(ii) + rVec(ii)*y_circle;
        plot(x,y,'Color',[0.7,0.7,0])        
        if labelcircles
            text(xCenVec(ii)+rVec(ii),yCenVec(ii),num2str(numbering(ii)),'FontSize',textSize,'Color',textColor);
        end
    end
    hold off;
end