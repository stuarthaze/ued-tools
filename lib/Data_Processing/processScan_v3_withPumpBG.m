function ScanData = processScan_v3_withPumpBG(dataDir,ROIs3d,IMG_BG_det,IMG_BG_pump, imgCrop,outlierThreshold)

imgSize = size(IMG_BG_det);
[NY_full,NX_full] = size(IMG_BG_det);
[NY_crop,NX_crop,Nspots] = size(ROIs3d);
backupDirectories = dir(fullfile(dataDir,'Backup*'));
% Read first directory to find number of frames for pre-allocation of
% memory
directoryName_1 = [dataDir, '\', backupDirectories(1).name];
filesPP_1 = dir(fullfile(directoryName_1, '*Pump+Probe*'));
numFrames = length(filesPP_1);
NT = length(backupDirectories);
nPixSpot = sum(sum(ROIs3d(:,:,1)));
% Pre-allocate
IMGs_P  = zeros(NY_full, NX_full, NT);
IMGs_PP = zeros(NY_full, NX_full, NT);
  
imgs_t_P  = zeros(NY_full,NX_full,numFrames);
imgs_t_PP = zeros(NY_full,NX_full,numFrames);

usedFrames_P = zeros(NT,numFrames);
usedFrames_PP = zeros(NT,numFrames);

BG_3d_det = zeros(NY_full,NX_full,numFrames);
BG_3d_pump = zeros(NY_full,NX_full,numFrames);

for iFrame = 1:numFrames
    BG_3d_det(:,:,iFrame) = IMG_BG_det;
    BG_3d_pump(:,:,iFrame) = IMG_BG_pump;
end

spotIntensities_P = zeros(NT,Nspots,numFrames);
spotIntensities_PP = zeros(NT,Nspots,numFrames);
Intensities_P = zeros(NT,Nspots);
Intensities_PP = zeros(NT,Nspots);
Uncertainties_P = zeros(NT,Nspots);
Uncertainties_PP = zeros(NT,Nspots);
imageSumsAll_P = zeros(NT,numFrames);
imageSumsAll_PP = zeros(NT,numFrames);

for tIndex = 1:NT
    fprintf(['Reading time point ',num2str(tIndex), ' of ', num2str(NT),'\n'])
    
    directoryName = [dataDir, '\', backupDirectories(tIndex).name];
    imgs_t_P = readImagesInFolder(directoryName,'* Probe*',imgSize)-BG_3d_det;
    imgs_t_PP = readImagesInFolder(directoryName,'*Pump+Probe*',imgSize)-BG_3d_pump;

    fprintf(['Processing time point ',num2str(tIndex), ' of ', num2str(NT),'\n'])
    
    files_P = dir(fullfile(directoryName, '* Probe*'));
    files_PP = dir(fullfile(directoryName, '*Pump+Probe*'));
    searchPattern = '\d*(?=fs)';
    timeDelayCell = regexp(files_PP(1).name, searchPattern, 'match');
    timeDelays(tIndex) = str2double(timeDelayCell{1});
    
    disp(['numFrames = ',num2str(numFrames)]);
    
    for frame = 1:numFrames
        for spotIndex = 1:Nspots
            spotIntensities_P(tIndex,spotIndex,frame) = ...
                sum(sum(imgs_t_P(imgCrop.y,imgCrop.x,frame).*ROIs3d(:,:,spotIndex)))/nPixSpot;
            spotIntensities_PP(tIndex,spotIndex,frame) = ...
                sum(sum(imgs_t_PP(imgCrop.y,imgCrop.x,frame).*ROIs3d(:,:,spotIndex)))/nPixSpot;
        end
    end
    
    imgSums_P = sum(sum(imgs_t_P,2),1);
    imgSums_PP = sum(sum(imgs_t_PP,2),1);
    imageSumsAll_P(tIndex,:) = imgSums_P;
    imageSumsAll_PP(tIndex,:) = imgSums_PP;

    med_P = median(imgSums_P);
    med_PP = median(imgSums_PP);
    IQR_sum_P = quantile(imgSums_P,0.75)-quantile(imgSums_P,0.25);
    IQR_sum_PP = quantile(imgSums_PP,0.75)-quantile(imgSums_PP,0.25);
    IQR = 0.5*(IQR_sum_P + IQR_sum_PP);
    
    if IQR == 0
        disp('Warning IQR = 0');
    end

    goodframes_P = (imgSums_P <= (med_P + outlierThreshold*IQR)) & ...
        (imgSums_P >= (med_P - outlierThreshold*IQR));
    goodframes_PP = (imgSums_PP <= (med_PP + outlierThreshold*IQR)) & ...
        (imgSums_PP >= (med_PP - outlierThreshold*IQR));

    numGoodFrames_P = sum(goodframes_P);
    numGoodFrames_PP = sum(goodframes_PP);
    
    IMGs_P(:,:,tIndex) = mean(imgs_t_P(:,:,goodframes_P),3);
    IMGs_PP(:,:,tIndex) = mean(imgs_t_PP(:,:,goodframes_PP),3);


    Intensities_P(tIndex,:) = mean(spotIntensities_P(tIndex,:,goodframes_P),3);
    Intensities_PP(tIndex,:) = mean(spotIntensities_PP(tIndex,:,goodframes_PP),3);
    Uncertainties_P(tIndex,:) = std(spotIntensities_P(tIndex,:,goodframes_P),0,3)/sqrt(numGoodFrames_P);
    Uncertainties_PP(tIndex,:) = std(spotIntensities_PP(tIndex,:,goodframes_PP),0,3)/sqrt(numGoodFrames_PP);
%         
%     for spotIndex = 1:Nspots
%         Intensities_P(tIndex,spotIndex) = mean(spotIntensities_P(tIndex,spotIndex,goodframes_P),3);
%         Intensities_PP(tIndex,spotIndex) = mean(spotIntensities_PP(tIndex,spotIndex,goodframes_PP),3);
%         Uncertainties_P(tIndex,spotIndex) = std(spotIntensities_P(tIndex,spotIndex,goodframes_P),0,3)/sqrt(numGoodFrames_P);
%         Uncertainties_PP(tIndex,spotIndex) = std(spotIntensities_PP(tIndex,spotIndex,goodframes_PP),0,3)/sqrt(numGoodFrames_PP);
%     end
    numRejected_P(tIndex) = numFrames - numGoodFrames_P;
    numRejected_PP(tIndex) = numFrames - numGoodFrames_PP;
    IQR_T(tIndex) = IQR;
    usedFrames_P(tIndex,:) = goodframes_P;
    usedFrames_PP(tIndex,:) = goodframes_PP;
end
[timeDelaysSorted, sortingIndicies] = sort(timeDelays);
%--------------------------------------------------------------------------
% Check for problems based on IQR values
medianIQR = median(IQR_T);
IQR_IQR = iqr(IQR_T);
STD_IQR = std(IQR_T);
%--------------------------------------------------------------------------
% Statistics for all frames
medianImgSum_P = median(reshape(imageSumsAll_P,[NT*numFrames,1]));
medianImgSum_PP = median(reshape(imageSumsAll_PP,[NT*numFrames,1]));
iqrImgSum_P = iqr(reshape(imageSumsAll_P,[NT*numFrames,1]));
iqrImgSum_PP = iqr(reshape(imageSumsAll_PP,[NT*numFrames,1]));
%----------------------------------
IQRerrorThreshold = mean([medianIQR, iqrImgSum_P, iqrImgSum_PP]) + 2*IQR_IQR;
tIndexRange = 1:NT;
IQRerrorDetected = IQR_T > IQRerrorThreshold;
problemTimePoints = tIndexRange(IQRerrorDetected);
%--------------------------------------------------------------------------
% Reprocess timepoints with anomalous IQR
for tIndex = problemTimePoints
    fprintf(['Reprocessing time point ',num2str(tIndex), ' of ', num2str(NT),'\n'])
    
    directoryName = [dataDir, '\', backupDirectories(tIndex).name];
    imgs_t_P = readImagesInFolder(directoryName,'* Probe*',imgSize)-BG_3d_det;
    imgs_t_PP = readImagesInFolder(directoryName,'*Pump+Probe*',imgSize)-BG_3d_det;
    
%     files_P = dir(fullfile(directoryName, '* Probe*'));
    files_PP = dir(fullfile(directoryName, '*Pump+Probe*'));
    searchPattern = '\d*(?=fs)';
    timeDelayCell = regexp(files_PP(1).name, searchPattern, 'match');
    timeDelays(tIndex) = str2double(timeDelayCell{1});
       
    for frame = 1:numFrames
        for spotIndex = 1:Nspots
            spotIntensities_P(tIndex,spotIndex,frame) = ...
                sum(sum(imgs_t_P(imgCrop.y,imgCrop.x,frame).*ROIs3d(:,:,spotIndex)))/nPixSpot;
            spotIntensities_PP(tIndex,spotIndex,frame) = ...
                sum(sum(imgs_t_PP(imgCrop.y,imgCrop.x,frame).*ROIs3d(:,:,spotIndex)))/nPixSpot;
        end
    end
    
    imgSums_P = sum(sum(imgs_t_P,2),1);
    imgSums_PP = sum(sum(imgs_t_PP,2),1);
    imageSumsAll_P(tIndex,:) = imgSums_P;
    imageSumsAll_PP(tIndex,:) = imgSums_PP;

    newTimeIndex = tIndexRange(sortingIndicies == tIndex);
    if newTimeIndex == 1
        nearestTimePoints = sortingIndicies([0,1]+newTimeIndex);
        nAveraged = 2;
    elseif newTimeIndex == NT
        nearestTimePoints = sortingIndicies([-1,0]+newTimeIndex);
        nAveraged = 2;
    else
        nearestTimePoints = sortingIndicies([-1,0,1]+newTimeIndex);
        nAveraged = 3;
    end

    % recalculate medians including neighbouring time points
    med_P = median(reshape(imageSumsAll_P(nearestTimePoints,:),[numFrames*nAveraged,1]));
    med_PP = median(reshape(imageSumsAll_PP(nearestTimePoints,:),[numFrames*nAveraged,1]));


    goodframes_P = (imgSums_P <= (med_P + IQRerrorThreshold)) & ...
        (imgSums_P >= (med_P - IQRerrorThreshold));
    goodframes_PP = (imgSums_PP <= (med_PP + IQRerrorThreshold)) & ...
        (imgSums_PP >= (med_PP - IQRerrorThreshold));

    numGoodFrames_P = sum(goodframes_P);
    numGoodFrames_PP = sum(goodframes_PP);
    
    IMGs_P(:,:,tIndex) = mean(imgs_t_P(:,:,goodframes_P),3);
    IMGs_PP(:,:,tIndex) = mean(imgs_t_PP(:,:,goodframes_PP),3);

    Intensities_P(tIndex,:) = mean(spotIntensities_P(tIndex,:,goodframes_P),3);
    Intensities_PP(tIndex,:) = mean(spotIntensities_PP(tIndex,:,goodframes_PP),3);
    Uncertainties_P(tIndex,:) = std(spotIntensities_P(tIndex,:,goodframes_P),0,3)/sqrt(numGoodFrames_P);
    Uncertainties_PP(tIndex,:) = std(spotIntensities_PP(tIndex,:,goodframes_PP),0,3)/sqrt(numGoodFrames_PP);

    numRejected_P(tIndex) = numFrames - numGoodFrames_P;
    numRejected_PP(tIndex) = numFrames - numGoodFrames_PP;
    usedFrames_P(tIndex,:) = goodframes_P;
    usedFrames_PP(tIndex,:) = goodframes_PP;
end

IMG_diffs = IMGs_PP - IMGs_P;
ScanData.IMG_diffs = IMG_diffs(:,:,sortingIndicies);
ScanData.timeDelays = timeDelaysSorted;
ScanData.IMGs_P = IMGs_P(:,:,sortingIndicies);
ScanData.IMGs_PP = IMGs_PP(:,:,sortingIndicies);

ScanData.Intensities_P = Intensities_P(sortingIndicies,:);
ScanData.Intensities_PP = Intensities_PP(sortingIndicies,:);
ScanData.Uncertainties_P = Uncertainties_P(sortingIndicies,:);
ScanData.Uncertainties_PP = Uncertainties_PP(sortingIndicies,:);
ScanData.numRejectedFrames_P = numRejected_P(sortingIndicies);
ScanData.numRejectedFrames_PP = numRejected_PP(sortingIndicies);
ScanData.IQRs = IQR_T(sortingIndicies);
ScanData.IQR_reprocessingThreshold = IQRerrorThreshold;
ScanData.IntensitiesAllFrames_P = spotIntensities_P(sortingIndicies,:,:);
ScanData.IntensitiesAllFrames_PP = spotIntensities_PP(sortingIndicies,:,:);
ScanData.imageSumsAllFrames_P = imageSumsAll_P(sortingIndicies,:);
ScanData.imageSumsAllFrames_PP = imageSumsAll_PP(sortingIndicies,:);
ScanData.usedFrames_P = usedFrames_P(sortingIndicies,:);
ScanData.usedFrames_PP = usedFrames_PP(sortingIndicies,:);
ScanData.reprocessedTimePoints = IQRerrorDetected(sortingIndicies);


function IMGs = readImagesInFolder(DirectoryPath,string2match,img_size)
files = dir(fullfile(DirectoryPath, ['*',string2match,'*']));
numFiles = length(files);
IMGs = zeros([img_size, numFiles]);
for indx = 1:numFiles
    IMGs(:,:,indx) = double(DataIOLibrary.DataIO.ReadSpe(fullfile(DirectoryPath,files(indx).name)));
end
end

end