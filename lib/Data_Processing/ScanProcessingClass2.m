classdef ScanProcessingClass2
    
    properties
        dataDir
        backupDirs
        IMGs_P
        IMGs_PP
        IMGs_dI
        BG_Pump
        BG_Det
        beamstop
        ROIstack
        cropping
        nx
        nx_full
        ny
        ny_full
        nt
        nframes
        nroi
        numPixROI
        timeDelays
        centreRaw
        centre
        
        % Processing Options
        useDetBG
        useDynamicDetBG
        usePumpBG
        useDynamicPumpBG

        filterOptionsAveraged
        filterBackupFrames
        filterOptionsBackupFrames
        filterBackupProcessed
        filterOptionsBackupProcessed
        filterBackupDiffs
        filterOptionsBackupDiffs

        useIntensityCalibration
        CalibrationMatCropped

        % Reduced data
        Intensities_P       % [nt x nroi]
        Intensities_PP      % [nt x nroi]
        Intensities_dI      % [nt x nroi]
        RelativeChanges
        
        % Radial Background
        BG_radial1d
        Intensities_radialBG  % [nt x nroi]
        quantileRadialBG
        usedRadialBG
        usedBeamstop
        
        IntensitiesAllFrames_P  % [nt x nroi x nframes]
        IntensitiesAllFrames_PP % [nt x nroi x nframes]
        IntensitiesAllFrames_dI
        
        % Statistics
        Uncertainties_P     % [nt x nroi]
        Uncertainties_PP    % [nt x nroi]
        Uncertainties_dI    % [nt x nroi]
        Uncertainties_rel
        SNRmatrix
        rmsSNR_time
        rmsSNR_peaks
        
        imageSumsAll_P
        imageSumsAll_PP
        imageSumsAll_dI
        
        IQR_T
        IQR_allFrames_P
        IQR_allFrames_PP
        thresholdImageSum
        thresholdMaxOnly
        thresholdDifferenceOnly
        usedFrames
        numRejectedFrames
        
        IQR_reprocessingThreshold
        reprocessedTimePoints
    end %properties
%--------------------------------------------------------------   
    methods
        %Constructor
        function obj = ScanProcessingClass2(dataDir,ImageClass)
            obj.dataDir = dataDir;
            obj.cropping = ImageClass.cropping;
            obj.centreRaw = ImageClass.centreRaw;
            obj.centre = ImageClass.centre;
            obj.ROIstack = logical(ImageClass.ROIstack);
            obj.nx = size(obj.ROIstack,2);
            obj.ny = size(obj.ROIstack,1);
            obj.nx_full = ImageClass.imgsizeRaw(2);
            obj.ny_full = ImageClass.imgsizeRaw(1);
            obj.nroi = size(obj.ROIstack,3);
            obj.numPixROI = reshape(sum(sum(obj.ROIstack,1),2),[1,obj.nroi]);
            obj.thresholdImageSum = 3;
            % Set defaults for options
%             obj.useMedianFilter = 0;
            obj.useDetBG = 0;
            obj.useDynamicDetBG = 0;
            obj.usePumpBG = 0;
            obj.useDynamicPumpBG = 0;
            obj.useIntensityCalibration = 0;
            if ~isempty(ImageClass.beamstop)
                obj.beamstop = ImageClass.beamstop;
            end
        end 
        %------------------------------------------
        function obj = readStaticDetectorBG(obj,BG)
            if size(BG,1)==obj.ny_full && size(BG,2)==obj.nx_full
                obj.BG_Det = BG(obj.cropping.y,obj.cropping.x);
            elseif size(BG,1)==obj.ny && size(BG,2)==obj.nx
                obj.BG_Det = BG;
            else
                disp('Error in BG image size');
                return
            end
            obj.useDetBG = 1;
        end
        
        function obj = readStaticPumpBG(obj,BG)
            if size(BG,1)==obj.ny_full && size(BG,2)==obj.nx_full
                obj.BG_Pump = BG(obj.cropping.y,obj.cropping.x);
            elseif size(BG,1)==obj.ny && size(BG,2)==obj.nx
                obj.BG_Pump = BG;
            else
                disp('Error in BG image size');
                return
            end
            obj.usePumpBG = 1;
        end
        
        function obj = ReadCalibrationMat(obj,ImageClass)
            obj.useIntensityCalibration = 1;
            obj.CalibrationMatCropped = ImageClass.calibrationIMG(obj.cropping.y,obj.cropping.x);
        end
        
%--------------------------------------------------------------------------
        function obj = readScanAveragedImages(obj,options)
            if options.readTimeResolved
                disp('Reading Time-resolved images');
                [obj.IMGs_dI, timeDelaysTR] = obj.readAveragedImagesSorted('*Time-Resolved*');
                obj.nt = length(timeDelaysTR);
                obj.timeDelays = timeDelaysTR;
            end
            if options.readPump
                disp('Reading pump images');
                [obj.BG_Pump, ~] = obj.readAveragedImagesSorted('*Pump *');
            end
            if options.readProbe
                disp('Reading probe images');
                [obj.IMGs_P, timeDelaysProbe] = obj.readAveragedImagesSorted('* Probe*');
                ntProbe = length(timeDelaysProbe);
            end
            if options.readPumpProbe
                disp('Reading pump+probe images');
                [obj.IMGs_PP, timeDelaysPP] = obj.readAveragedImagesSorted('*Pump+Probe*');
                ntPP = length(timeDelaysPP);
                if isempty(obj.nt)
                    obj.nt = ntPP;
                    obj.timeDelays = timeDelaysPP;
                end
            end
            
            % Background subtractions
            if options.subtractDetectorBG
                for T = 1:ntProbe
                    obj.IMGs_P(:,:,T) = obj.IMGs_P(:,:,T) - obj.BG_Det;
                end
            end
            if options.subtractPumpBG
                if options.readPump
                    obj.IMGs_PP = obj.IMGs_PP - obj.BG_Pump;
                else
                    for T = 1:ntPP
                        obj.IMGs_PP(:,:,T) = obj.IMGs_PP(:,:,T) - obj.BG_Pump;
                    end
                end
            elseif options.subtractDetectorBG && ~isempty(obj.IMGs_PP)
                for T = 1:obj.nt
                    obj.IMGs_PP(:,:,T) = obj.IMGs_PP(:,:,T) - obj.BG_Det;
                end
            end
            
            % Calibration
            if obj.useIntensityCalibration
                if isempty(obj.CalibrationMatCropped)
                    disp('Warning: no calibration data');
                elseif options.readPumpProbe && options.readProbe
                    for T = 1:obj.nt
                        obj.IMGs_PP(:,:,T) = obj.IMGs_PP(:,:,T).*obj.CalibrationMatCropped;
                        obj.IMGs_P(:,:,T) = obj.IMGs_P(:,:,T).*obj.CalibrationMatCropped;
                    end
                end
            end
            
            if ~options.readTimeResolved && options.readPumpProbe && options.readProbe
                obj.IMGs_dI = obj.IMGs_PP - obj.IMGs_P;
            end
        end %readScanAveragedImages
        %------------------------------------------------------------------
        % Normalize
        function obj = normalizeAverages(obj)
            if isempty(obj.IMGs_P)
                disp('Error: no Probe images');
                return
            end
            if isempty(obj.IMGs_PP)
                disp('Error: no Pump+Probe images');
                return
            end
            if ~isempty(obj.beamstop)
                for iT = 1:obj.nt
                    obj.IMGs_P(:,:,iT) = (1-obj.beamstop).*obj.IMGs_P(:,:,iT);
                    obj.IMGs_PP(:,:,iT) = (1-obj.beamstop).*obj.IMGs_PP(:,:,iT);
                end
            end
            Psums = sum(sum(obj.IMGs_P,2),1);
            PPsums = sum(sum(obj.IMGs_PP,2),1);
            meanPsum = mean(Psums);
            meanPPsum = mean(PPsums);
            normP = meanPsum/Psums;
            normPP = meanPPsum/PPsums;
            for iT = 1:obj.nt
                obj.IMGs_P(:,:,iT) = obj.IMGs_P(:,:,iT)*normP(iT);
                obj.IMGs_PP(:,:,iT) = obj.IMGs_PP(:,:,iT)*normPP(iT);
                if ~isempty(obj.Intensities_P) && ~isempty(obj.Intensities_PP)
                    obj.Intensities_P(iT,:) = obj.Intensities_P(iT,:)*normP(iT);
                    obj.Intensities_PP(iT,:) = obj.Intensities_PP(iT,:)*normPP(iT);
                end
            end
            obj.IMGs_dI = obj.IMGs_PP - obj.IMGs_P;
            if ~isempty(obj.Intensities_P) && ~isempty(obj.Intensities_PP)
                obj.Intensities_dI = obj.Intensities_PP - obj.Intensities_P;
                obj = obj.calcRelativeChanges();
            end
        end
        %------------------------------------------------------------------
        function obj = extractROIvalsFromImages(obj)
            if ~isempty(obj.IMGs_dI)
                obj.Intensities_dI = obj.imgStack2roivals(obj.IMGs_dI);
            end
            if ~isempty(obj.IMGs_P)
                obj.Intensities_P = obj.imgStack2roivals(obj.IMGs_P);
            end
            if ~isempty(obj.IMGs_PP)
                obj.Intensities_PP = obj.imgStack2roivals(obj.IMGs_PP);
            end
            if ~isempty(obj.Intensities_P) && ~isempty(obj.Intensities_dI)
                obj = obj.calcRelativeChanges();
            end
        end
        
        function ROIvals = img2roivals(obj,IMG)
            ROIvals = zeros(1,obj.nroi);
            for ii = 1:obj.nroi
                ROIvals(ii) = sum(IMG(obj.ROIstack(:,:,ii)))/obj.numPixROI(ii);
            end
        end
        function Vals = imgStack2roivals(obj,IMGstack)
            nImages = size(IMGstack,3);
            Vals = zeros([nImages,obj.nroi]);
            for indxImg = 1:nImages
                Vals(indxImg,:) = obj.img2roivals(IMGstack(:,:,indxImg));
            end
        end
          
        %------------------------------------------------------------------
        function obj = processScanFromBackups(obj)
            obj.backupDirs = dir(fullfile(obj.dataDir,'Backup*'));
            % Read first directory to find number of frames for pre-allocation of memory
            directoryName_1 = [obj.dataDir, '\', obj.backupDirs(1).name];
            filesPP_1 = dir(fullfile(directoryName_1, '*Pump+Probe*'));
            obj.nframes = length(filesPP_1);
            NT = length(obj.backupDirs);
            obj.nt = NT;
            % Pre-allocate - Averaged images
            obj.IMGs_P  = zeros(obj.ny, obj.nx, NT);
            obj.IMGs_PP = zeros(obj.ny, obj.nx, NT);
            if obj.useDynamicPumpBG
                obj.BG_Pump = zeros(obj.ny, obj.nx, NT);
            end

            % Preallocate          
            obj.usedFrames = zeros(NT,obj.nframes);
            % Preallocate memory for individual frames
            framesDet = zeros(obj.ny,obj.nx,obj.nframes);
            framesPump = zeros(obj.ny,obj.nx,obj.nframes);
%             framesProbe = zeros(obj.ny_full,obj.nx_full,obj.nframes);
%             framesPP = zeros(obj.ny_full,obj.nx_full,obj.nframes);
            
            % Set Static Detector BG
            if obj.useDetBG
                if ~obj.useDynamicDetBG
                    if isempty(obj.BG_Det)
                        disp('Load detector background, BG_Det');
                        return
                    else
                        BG = obj.BG_Det;
%                         if obj.useMedianFilter
%                             BG = medfilt2(BG,obj.medianFilterWindow(1:2));
%                         end
                        if obj.filterBackupFrames
                            BG = outlierFilterGeneral(BG,obj.filterOptionsBackupFrames);
                        end

                        for iframe = 1:obj.nframes
                            framesDet(:,:,iframe) = BG;
                        end
                    end
                end
            else
                disp('Warning - no detector BG subtraction')
            end
            % Set static Pump BG
            if obj.usePumpBG
                if ~obj.useDynamicPumpBG
                    if isempty(obj.BG_Pump)
                        disp("Load BG_pump");
                        return
                    else
                        BG = obj.BG_Pump;
%                         if obj.useMedianFilter
%                             BG = medfilt2(BG,obj.medianFilterWindow(1:2));
%                         end

                        if obj.filterBackupFrames
                            BG = outlierFilterGeneral(BG,obj.filterOptionsBackupFrames);
                        end
                        for iframe = 1:obj.nframes
                            framesPump(:,:,iframe) = BG;
                        end
                    end
                end
            else
                framesPump = framesDet;
            end
            
            % Preallocate
            obj.IntensitiesAllFrames_P = zeros(NT,obj.nroi,obj.nframes);
            obj.IntensitiesAllFrames_PP = zeros(NT,obj.nroi,obj.nframes); 
            obj.IntensitiesAllFrames_dI = zeros(NT,obj.nroi,obj.nframes); 
            obj.Intensities_P = zeros(NT,obj.nroi);
            obj.Intensities_PP = zeros(NT,obj.nroi);
            obj.Intensities_dI = zeros(NT,obj.nroi);
            obj.Uncertainties_P = zeros(NT,obj.nroi);
            obj.Uncertainties_PP = zeros(NT,obj.nroi);
            obj.Uncertainties_dI = zeros(NT,obj.nroi);
            obj.imageSumsAll_P = zeros(NT,obj.nframes);
            obj.imageSumsAll_PP = zeros(NT,obj.nframes);
            obj.imageSumsAll_dI = zeros(NT,obj.nframes);

            % ------------- Loop over timepoints ----------------
            for tIndex = 1:NT
                fprintf(['Reading time point ',num2str(tIndex), ' of ', num2str(NT),'\n'])
                directoryName = [obj.dataDir, '\', obj.backupDirs(tIndex).name];
                % Read all frames at tIndex
                disp('Reading Probe');
                framesProbe = obj.readBackupFramesSorted(directoryName,'* Probe*');
                disp('Reading Pump+Probe');
                framesPP = obj.readBackupFramesSorted(directoryName,'*Pump+Probe*');
                if obj.useDynamicPumpBG
                    disp('Reading Pump');
                    framesPump = obj.readBackupFramesSorted(directoryName,'*Pump *');
                end
                if obj.useDynamicDetBG
                    disp('Warning - dynamic det BG not tested');
                    framesDet = obj.readBackupFramesSorted(directoryName,'* Detector*');
                end

                % Display status
                fprintf(['Processing time point ',num2str(tIndex), ' of ', num2str(NT),'\n'])
% 
%                 % Filter
%                 if obj.useMedianFilter
%                     framesProbe = medfilt3(framesProbe,obj.medianFilterWindow);
%                     framesPP = medfilt3(framesPP,obj.medianFilterWindow);
%                     if obj.useDynamicDetBG
%                         framesDet = medfilt3(framesDet,obj.medianFilterWindow);
%                     end
%                     if obj.useDynamicPumpBG
%                         framesPump = medfilt3(framesPump,obj.medianFilterWindow);
%                     end
%                 end
%                 
%                 if obj.useOutlierFilter
%                     framesProbe = outlierFilter3d(framesProbe,obj.outlierThreshold);
%                     framesPP = outlierFilter3d(framesPP,obj.outlierThreshold);
%                     framesDet = outlierFilter3d(framesDet,obj.outlierThreshold);
%                     framesPump = outlierFilter3d(framesPump,obj.outlierThreshold);
%                 end

                if obj.filterBackupFrames
                    framesProbe = outlierFilterGeneral(framesProbe,obj.filterOptionsBackupFrames);
                    framesPP = outlierFilterGeneral(framesPP,obj.filterOptionsBackupFrames);
                    if obj.useDynamicPumpBG
                        framesPump = outlierFilterGeneral(framesPump,obj.filterOptionsBackupFrames);
                    end
                    if obj.useDynamicDetBG
                        framesDet = outlierFilterGeneral(framesDet,obj.filterOptionsBackupFrames);
                    end
                end

                % Subtract BGs
                imgsProbe = framesProbe - framesDet;
                imgsPP = framesPP - framesPump;

                % Calibration
                if obj.useIntensityCalibration
                    for frame = 1:obj.nframes
                        imgsProbe(:,:,frame) = imgsProbe(:,:,frame).*obj.CalibrationMatCropped;
                        imgsPP(:,:,frame) = imgsPP(:,:,frame).*obj.CalibrationMatCropped;
                    end
                end
                
                imgsDiffs = imgsPP - imgsProbe;
                if obj.filterBackupDiffs
                    for iter = 1:obj.filterBackupDiffs
                        imgsDiffs = outlierFilterGeneral(imgsDiffs,obj.filterOptionsBackupDiffs);
                    end
                end
                
                if obj.filterBackupProcessed
                    for iter = 1:obj.filterBackupProcessed
                        imgsProbe = outlierFilterGeneral(imgsProbe,obj.filterOptionsBackupProcessed);
                        imgsPP = outlierFilterGeneral(imgsPP,obj.filterOptionsBackupProcessed);
                    end
                end


                % Read timedelay
                files_PP = dir(fullfile(directoryName, '*Pump+Probe*'));
                searchPattern = '\d*(?=fs)';
                timeDelayCell = regexp(files_PP(1).name, searchPattern, 'match');
                obj.timeDelays(tIndex) = str2double(timeDelayCell{1});

                % Extract intensities from ROIs
                for frame = 1:obj.nframes
                    obj.IntensitiesAllFrames_P(tIndex,:,frame) = obj.img2roivals(imgsProbe(:,:,frame));
                    obj.IntensitiesAllFrames_PP(tIndex,:,frame) = obj.img2roivals(imgsPP(:,:,frame));
                    obj.IntensitiesAllFrames_dI(tIndex,:,frame) = obj.img2roivals(imgsDiffs(:,:,frame));
                end
                
                % Check for errors based on roi vals
                errorThreshold = 3;
                ratioP = zeros(obj.nroi,obj.nframes);
                ratioPP = zeros(obj.nroi,obj.nframes);
                errorsROIsP = zeros(obj.nroi,obj.nframes);
                errorsROIsPP = zeros(obj.nroi,obj.nframes);
                
                for roi = 1:obj.nroi
                    roiP = obj.IntensitiesAllFrames_P(tIndex,roi,:);
                    medianP = median(roiP);
                    ratioP(roi,:) = roiP/medianP;
                    roiPP = obj.IntensitiesAllFrames_PP(tIndex,roi,:);
                    medianPP = median(roiPP);
                    ratioPP(roi,:) = roiPP/medianPP;
                    %find errors
                    iqrP = quantile(ratioP(roi,:)-1,0.75) - quantile(ratioP(roi,:)-1,0.25);
                    errorsROIsP = abs(ratioP(roi,:)-1) > iqrP*errorThreshold;
                    
                    iqrPP = quantile(ratioPP(roi,:)-1,0.75) - quantile(ratioPP(roi,:)-1,0.25);
                    errorsROIsPP = abs(ratioPP(roi,:)-1) > iqrPP*errorThreshold;
                end

% old version                errorsROIsPP = abs(ratioPP-1) >= errorThreshold;
                errorsROIs = errorsROIsP | errorsROIsPP;
                goodframesROIs = ~logical(sum(errorsROIs,1));

                % Check for errors based on sum(img)
                
                imgSums_P  = sum(imgsProbe,[1,2]); 
                imgSums_PP = sum(imgsPP,[1,2]);
                imgSums_dI = sum(imgsDiffs,[1,2]);
                
                obj.imageSumsAll_P(tIndex,:)  = imgSums_P; 
                obj.imageSumsAll_PP(tIndex,:) = imgSums_PP;
                obj.imageSumsAll_dI(tIndex,:) = imgSums_dI;

                % Median values at current time point
                med_P = median(obj.imageSumsAll_P(tIndex,:));
                med_PP = median(obj.imageSumsAll_PP(tIndex,:));
                med_diff = median(obj.imageSumsAll_dI(tIndex,:));
                
                % IQR over all available time points
                IQR_sum_P  = iqr(obj.imageSumsAll_P(1:tIndex,:),'all');
                IQR_sum_PP = iqr(obj.imageSumsAll_PP(1:tIndex,:),'all');
                IQR_diff   = iqr(obj.imageSumsAll_dI(1:tIndex,:),'all');
                
% old version
%                 IQR_sum_P  = quantile(imgSums_P,0.75) -  quantile(imgSums_P,0.25);
%                 IQR_sum_PP = quantile(imgSums_PP,0.75) - quantile(imgSums_PP,0.25);
%                 IQR_diff = quantile(imgSums_dI,0.75) - quantile(imgSums_dI,0.25);
                
                IQR = (IQR_sum_P + IQR_sum_PP + IQR_diff)/3;
                obj.IQR_T(tIndex) = IQR;

                if IQR == 0 || IQR_diff == 0 || IQR_sum_P==0 || IQR_sum_PP ==0
                    disp('Warning IQR = 0');
                    IQR
                    IQR_sum_P
                    IQR_sum_PP
                    IQR_diff
                    IQR = Inf;
                end

                if obj.thresholdMaxOnly
                    goodframes_P = (imgSums_P <= (med_P + obj.thresholdImageSum*IQR_sum_P));
                    goodframes_PP = (imgSums_PP <= (med_PP + obj.thresholdImageSum*IQR_sum_PP));
                    goodframes = goodframes_P & goodframes_PP;
                elseif obj.thresholdDifferenceOnly
                    goodframes_diff = (imgSums_dI <= (med_diff + obj.thresholdImageSum*IQR_diff)) & ...
                                      (imgSums_dI >= (med_diff - obj.thresholdImageSum*IQR_diff));
                    goodframes = goodframes_diff;
                    IQR = IQR_diff;
                else
                    goodframes_P = (imgSums_P <= (med_P + obj.thresholdImageSum*IQR_sum_P)) & ...
                                   (imgSums_P >= (med_P - obj.thresholdImageSum*IQR_sum_P));
                    goodframes_PP = (imgSums_PP <= (med_PP + obj.thresholdImageSum*IQR_sum_PP)) & ...
                                    (imgSums_PP >= (med_PP - obj.thresholdImageSum*IQR_sum_PP));
                    goodframes_diff = (imgSums_dI <= (med_diff + obj.thresholdImageSum*IQR_diff)) & ...
                                      (imgSums_dI >= (med_diff - obj.thresholdImageSum*IQR_diff));
                    goodframes = goodframes_P & goodframes_PP & goodframes_diff;
                end
                goodframes = reshape(goodframes,[1,obj.nframes]);
                goodframes = goodframes & goodframesROIs;

                numGoodFrames = sum(goodframes);
                
                % Store results in properties
                obj.IMGs_P(:,:,tIndex)   = mean(imgsProbe(:,:,goodframes),3);
                obj.IMGs_PP(:,:,tIndex)  = mean(imgsPP(:,:,goodframes),3);
                obj.IMGs_dI(:,:,tIndex)  = mean(imgsDiffs(:,:,goodframes),3);
                obj.Intensities_P(tIndex,:) = mean(obj.IntensitiesAllFrames_P(tIndex,:,goodframes),3);
                obj.Intensities_PP(tIndex,:)= mean(obj.IntensitiesAllFrames_PP(tIndex,:,goodframes),3);
                obj.Intensities_dI(tIndex,:)= mean(obj.IntensitiesAllFrames_dI(tIndex,:,goodframes),3);
                obj.Uncertainties_P(tIndex,:)  = std(obj.IntensitiesAllFrames_P(tIndex,:,goodframes),0,3)/sqrt(numGoodFrames);
                obj.Uncertainties_PP(tIndex,:) = std(obj.IntensitiesAllFrames_PP(tIndex,:,goodframes),0,3)/sqrt(numGoodFrames);
                obj.Uncertainties_dI(tIndex,:) = std(obj.IntensitiesAllFrames_dI(tIndex,:,goodframes),0,3)/sqrt(numGoodFrames);
                obj.numRejectedFrames(tIndex) = obj.nframes - numGoodFrames;
                
                obj.usedFrames(tIndex,:) = goodframes;
                if obj.useDynamicPumpBG
                    obj.BG_Pump(:,:,tIndex) = mean(framesPump(:,:,goodframes),3);
                end

            end %loop over T
            obj = obj.sortDataByTime();
            obj = obj.calcRelativeChanges();
        end %processScanFromBackups
        %------------------------------------------------------------------
        
        function obj = sortDataByTime(obj)
            [obj.timeDelays, sortingIndices] = sort(obj.timeDelays);
            obj.IMGs_P  = obj.IMGs_P(:,:,sortingIndices);
            obj.IMGs_PP = obj.IMGs_PP(:,:,sortingIndices);
            obj.IMGs_dI = obj.IMGs_dI(:,:,sortingIndices);
            obj.Intensities_P  = obj.Intensities_P(sortingIndices,:);
            obj.Intensities_PP = obj.Intensities_PP(sortingIndices,:);
            obj.Intensities_dI = obj.Intensities_dI(sortingIndices,:);
            obj.Uncertainties_P  = obj.Uncertainties_P(sortingIndices,:);
            obj.Uncertainties_PP = obj.Uncertainties_PP(sortingIndices,:);
            obj.Uncertainties_dI = obj.Uncertainties_dI(sortingIndices,:);
            obj.numRejectedFrames = obj.numRejectedFrames(sortingIndices);
            obj.IQR_T = obj.IQR_T(sortingIndices);
            obj.usedFrames = obj.usedFrames(sortingIndices,:);
            obj.imageSumsAll_P = obj.imageSumsAll_P(sortingIndices,:);
            obj.imageSumsAll_PP = obj.imageSumsAll_PP(sortingIndices,:);
            obj.imageSumsAll_dI = obj.imageSumsAll_dI(sortingIndices,:);
            if obj.useDynamicPumpBG
                obj.BG_Pump = obj.BG_Pump(:,:,sortingIndices);
            end
        end
        
        %------------------------------------------------------------------
        % Radial Background
        function obj = calcRadialBGs(obj,quantile)
            if nargin == 2
                obj.quantileRadialBG = quantile;
            elseif isempty(obj.quantileRadialBG)
                disp('Error: Quantile value is missing');
            end
            obj.Intensities_radialBG = zeros([obj.nt,obj.nroi]);
            if ~isempty(obj.beamstop)
                block = obj.beamstop;
                obj.usedBeamstop = 1;
            else
                block = zeros([obj.ny,obj.nx]);
                obj.usedBeamstop = 0;
            end
            for T = 1:obj.nt
                disp(['Calculating radial BG no. ',num2str(T)]);
                [BGrad, obj.BG_radial1d(T,:)] = radialQuantile_BG(obj.IMGs_P(:,:,T),obj.centre,block,obj.quantileRadialBG);
                obj.Intensities_radialBG(T,:) = obj.img2roivals(BGrad);
            end
        end
     
        
        %------------------------------------------------------------------
        % Statistics
        %------------------------------------------------------------------
        function obj = calcRelativeChanges(obj,subtractRadialBG)
            I_probe = obj.Intensities_P;
            if nargin == 2 && subtractRadialBG && ~isempty(obj.Intensities_radialBG)
                I_probe = obj.Intensities_P - obj.Intensities_radialBG;
            end
            if min(min(I_probe)) < 0
                disp('Error in BG subtraction: I < 0');
            end
            obj.RelativeChanges = obj.Intensities_dI./I_probe;
            if ~isempty(obj.Uncertainties_dI)
                obj.Uncertainties_rel = obj.Uncertainties_dI./I_probe;
            end
        end
                
        function obj = calculateScanStats(obj)
            obj.IQR_allFrames_P = iqr(reshape(obj.imageSumsAll_P,[obj.nt*obj.nframes,1]));
            obj.IQR_allFrames_PP = iqr(reshape(obj.imageSumsAll_PP,[obj.nt*obj.nframes,1]));
            obj.SNRmatrix = obj.Intensities_dI./obj.Uncertainties_dI;
            obj.rmsSNR_time = rms(obj.SNRmatrix,2);
            obj.rmsSNR_peaks = rms(obj.SNRmatrix,1);
        end
        
        %----------------------------------------------------
        % Statistics for sample decay
        function showPeaksAndBG(obj,nPeaks,weakStrongCutoff)
            I_probe = obj.Intensities_P - obj.Intensities_radialBG;
            Istrong = mean(I_probe(:,1:weakStrongCutoff),2);
            Iweak = mean(I_probe(:,(weakStrongCutoff+1):nPeaks),2);
            Ibg = mean(obj.Intensities_radialBG,2);
            
            Istrong = Istrong./Istrong(1);
            Iweak = Iweak./Iweak(1);
            Ibg = Ibg./Ibg(1);
            
            figure(); plot([Istrong,Iweak,Ibg]);
            legend('Strong','Weak','BG');
            ylabel('Relative mean intensity');
            xlabel('Time index');
        end
        
        %------------------------------------------------------------------
        % Filtering Methods
        %------------------------------------------------------------------
        function obj = filterImagesMedian(obj,MxNxT)
            if ~isempty(obj.IMGs_dI)
                obj.IMGs_dI = medfilt3(obj.IMGs_dI,MxNxT);
            end
            if ~isempty(obj.IMGs_P)
                obj.IMGs_P = medfilt3(obj.IMGs_P,MxNxT);
            end
            if ~isempty(obj.IMGs_PP)
                obj.IMGs_PP = medfilt3(obj.IMGs_PP,MxNxT);
            end
        end
%         
%         function obj = setFilterOptionsAveraged(obj,options)
%             obj.filterOptionsAveraged = options;
%         end
        
        function obj = setFilterOptionsBackups(obj,options)
            obj.filterOptionsBackups = options;
        end
        
        function obj = filterImagesOutliers(obj,options)
            obj.filterOptionsAveraged = options;

            if ~isempty(obj.IMGs_dI)
                obj.IMGs_dI = outlierFilterGeneral(obj.IMGs_dI,options);
            end
            if ~isempty(obj.IMGs_P)
                obj.IMGs_P = outlierFilterGeneral(obj.IMGs_P,options);
            end
            if ~isempty(obj.IMGs_PP)
                obj.IMGs_PP = outlierFilterGeneral(obj.IMGs_PP,options);
            end
        end
  
        %------------------------------------------------------------------
        % Display results
        %------------------------------------------------------------------
        function showUsedFrames(obj)
            figure(); imagesc(obj.usedFrames);
            title('Frames Used');
            xlabel('Frame number'); ylabel('Time index');
        end
        
        function showDifferenceImages(obj,options)
            SCALE = [-options.scale,options.scale];
            nRows = options.nRows;
            nCols = options.nCols;
            if options.cropimg
                xcrop = options.xcrop;
                ycrop = options.ycrop;
            else
                xcrop = 1:obj.nx;
                ycrop = 1:obj.ny;
            end
            numberFigures =ceil(obj.nt/(nRows*nCols));
            imgsPerFig = nRows*nCols;
            for fig = 1:numberFigures
                figure();
                for row = 1:nRows
                    for col = 1:nCols
                        iImg = (fig-1)*imgsPerFig + (row-1)*nCols + col;
                        if iImg <= obj.nt
                            IMG = obj.IMGs_dI(ycrop,xcrop,iImg);
                            subplot(nRows,nCols,(row-1)*nCols + col);
                            if options.medianFilter
                                n = options.medianFilter;
                                IMG = medfilt2(IMG,[n,n]);
                            end
                            if options.movingAverageFilter
                                IMG = fMovingAverageFilter(IMG,obj.movingAverageFilter);
                            end
                            if options.gaussFilter
                                IMG = imgaussfilt(IMG,options.gaussFilter);
                            end
                            imagesc(IMG,SCALE);
                            pbaspect([length(xcrop),length(ycrop), 1])
                            colormap('gray');
                            title(num2str(obj.timeDelays(iImg)));
                        else
                            return
                        end
                    end
                end
            end            
        end %showDifferenceImages
        
        function IMGout = showAverageDifferenceImage(obj,imgIndices,SCALE)
            nIMGs = length(imgIndices);
            img = mean(obj.IMGs_dI(:,:,imgIndices),3);
            figure();
            imagesc(img,SCALE);
            pbaspect([obj.nx,obj.ny,1]); colormap('gray');
            if nIMGs == 1
                title(['img ',num2str(imgIndices),',  t = ',num2str(obj.timeDelays(imgIndices))]);
            else
                title(['Average of ',num2str(nIMGs),' images: ',num2str(imgIndices(1)),' - ',num2str(imgIndices(end))]);
            end
            IMGout = img;
        end
        %------------------------------------------------------------------------

        function showIntensityChanges(obj)
            showTRdata(obj.Intensities_dI);
            title('ROI intensity changes');
        end
        
        function showRelativeChanges(obj)
            showTRdata(obj.RelativeChanges);
            title('Relative Changes');
        end
        
        function showMeanIntensityChangesVsTime(obj,rois)
            if (nargin < 2) || isempty(rois)
                rois = 1:obj.nroi;
            end             
            dIvsT = mean(abs(obj.Intensities_dI(:,rois)),2);
            figure(); plot(obj.timeDelays*1e-3,dIvsT);
            xlabel('Time / ps'); title('Mean ROI intensity changes');
        end

        function showSNR_time(obj)
            figure();
            plot(obj.timeDelays*1e-3, obj.rmsSNR_time); 
            xlabel('Time / ps'); title('RMS SNR');
        end
            
        function showSNRmatrix(obj)
            showTRdata(obj.SNRmatrix);
            title('SNR matrix');
        end
        
        function showImageSumStatistics(obj)
            figure();
            plot(1:obj.nt, obj.IQR_T, 'k');
            title('IQR values of image sums'); xlabel('Time index');
            if ~isempty(obj.IQR_reprocessingThreshold)
                hline = refline(0,obj.reprocessingThreshold);
                text(1,1.1*obj.IQR_reprocessingThreshold, 'Reprocessing Threshold');
            end
            if ~isempty(obj.IQR_allFrames_P) && ~isempty(obj.IQR_allFrames_PP)
                hline = refline(0,obj.IQR_allFrames_P); hline.Color = 'b';
                t = text(0,0,'- Probe global IQR'); t.Color = 'b'; t.Units = 'normalized'; t.Position = [0.01, 0.95];
                hline = refline(0,obj.IQR_allFrames_PP); hline.Color = 'r';
                t = text(0,0,'- Pump+Probe global IQR'); t.Color = 'r'; t.Units = 'normalized'; t.Position = [0.01, 0.90];
            end
        end
        
        %----- Show ROI vs time --------
        function showRelativeChanges_errorbar(obj,rois)
            for ii = 1:length(rois)
                figure();
                roi = rois(ii);
                errorbar(obj.timeDelays,obj.RelativeChanges(:,roi),obj.Uncertainties_rel(:,roi));
                refline(0,0);
                title(['ROI #',num2str(roi)]);
            end
        end
        
        function saveRelativeChanges_errorbar(obj,rois)
            newDirectoryName = 'ROI Relative Changes Errorbar';
            if ~isdir(['.\',newDirectoryName])
                mkdir(newDirectoryName);
            end
            figure();
            for ii = 1:length(rois)
                roi = rois(ii);
                errorbar(obj.timeDelays,obj.RelativeChanges(:,roi),obj.Uncertainties_rel(:,roi));
                refline(0,0);
                title(['ROI #',num2str(roi)]);
                filename2save = ['plot_ROI_',num2str(roi),'.jpg'];
                print(filename2save,'-djpeg');
                movefile(filename2save,newDirectoryName); 
            end
        end

        
        %------------------------------------------------------------------
        % Read data
        function IMGs = readBackupFramesSorted(obj,DirectoryPath,string2match)
            files = dir(fullfile(DirectoryPath,string2match));
            numFiles = length(files);
            IMGs = zeros([obj.ny, obj.nx, numFiles]);
            for indx = 1:numFiles
                frameNumber = str2double(regexp(files(indx).name,'(?<=frame )\d*','match'));
                fullframe = double(DataIOLibrary.DataIO.ReadSpe(fullfile(DirectoryPath,files(indx).name)));
                IMGs(:,:,frameNumber) = fullframe(obj.cropping.y,obj.cropping.x); 
            end
        end
        
        function [IMGs, timeDelaysSorted] = readAveragedImagesSorted(obj,string2match)
            files = dir(fullfile(obj.dataDir, string2match));
            numFiles = length(files);
            IMGs = zeros([obj.ny,obj.nx,numFiles]);
            timedelays = zeros([1,numFiles]);
            for indx = 1:numFiles
                disp(['Image number ',num2str(indx),' of ',num2str(numFiles)]);
                timedelays(indx) = str2double(regexp(files(indx).name,'\d*(?=fs)','match'));
                fullframe = double(DataIOLibrary.DataIO.ReadSpe(fullfile(obj.dataDir,files(indx).name)));
                IMGs(:,:,indx) = fullframe(obj.cropping.y,obj.cropping.x);
            end
            [timeDelaysSorted, sortingIndices] = sort(timedelays);
            IMGs = IMGs(:,:,sortingIndices);
            if obj.nframes > 1
                IMGs = IMGs/obj.nframes;
            end
        end
    end %methods  
end
%-------------------------------
%      Utility Functions
%-------------------------------

function showTRdata(DATA)
    load('MyColormap_BWR','Colormap_BlueWhiteRed'); 
    ndat = size(DATA,1)*size(DATA,2);
    scale = quantile(reshape(abs(DATA),[ndat,1]),(1-2/ndat));
    SCALE = [-scale,scale];
    figname = figure(); imagesc(DATA,SCALE);
    set(figname,'Colormap',Colormap_BlueWhiteRed);
    xlabel('Peak Index'); ylabel('Time Index');
    colorbar();
end
% 
% function [IMGout,outliers] = outlierFilter(IMGin,threshold)
% % Replaces pixels in IMGin with median value of neigbours if it differs by 
% % more than threshold*std2(IMG-IMGmed)
%     IMGmed = medfilt2(IMGin);
%     IMGdiff = IMGin-IMGmed;
%     stdIMG = std2(IMGdiff);
%     outliers = (IMGdiff >= threshold*stdIMG) | (IMGdiff <= -threshold*stdIMG);
%     IMGout = IMGin;
%     IMGout(outliers) = IMGmed(outliers);
% end
% 
% function [IMGsOUT,outliers3d] = outlierFilter3d(IMGs3d,threshold)
%     [ny,nx,nz] = size(IMGs3d);
%     if nz == 1
%         [IMGsOUT,outliers3d] = outlierFilter(IMGs3d,threshold);
%     else
%         IMGsOUT = IMGs3d;
%         outliers3d = zeros([ny,nx,nz],'logical');
%         for z = 1:nz
%             [IMGsOUT(:,:,z),outliers3d(:,:,z)] = outlierFilter(IMGs3d(:,:,z),threshold);
%         end
%     end
% end