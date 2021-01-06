function scanData = processScan_TimeResolvedImages(DirectoryPath,imgCrop,ROIs_3d)
dataIOAssembly = NET.addAssembly('\\win.desy.de\home\hayess\My Documents\MATLAB\lib\Data_Processing\DataIOLib.dll');
files = dir(fullfile(DirectoryPath, '*Time-Resolved*'));
numFiles = length(files);
nxCrop = length(imgCrop.x);
nyCrop = length(imgCrop.y);
nROIs = size(ROIs_3d,3);
npixROI = reshape(sum(sum(ROIs_3d,2),1),[1,nROIs]);
IMGstack_cropped = zeros(nyCrop,nxCrop,numFiles);
dI = zeros(numFiles,nROIs);
for indx = 1:numFiles
    disp(['processing image ',num2str(indx)]);
    timeDelays(indx) = str2double(regexp(files(indx).name, '\d*(?=fs)', 'match'));
    IMG = double(DataIOLibrary.DataIO.ReadSpe(fullfile(DirectoryPath,files(indx).name)));
    IMGstack_cropped(:,:,indx) = IMG(imgCrop.y,imgCrop.x);
    for roi = 1:nROIs
        dI(indx,roi) = sum(sum(ROIs_3d(:,:,roi).*IMG(imgCrop.y,imgCrop.x)))/npixROI(roi);
    end
end

[timeDelaysSorted, sortingIndicies] = sort(timeDelays);
scanData.timeDelays = timeDelaysSorted;
scanData.dI_ROIs = dI(sortingIndicies,:);
scanData.IMGs_TR = IMGstack_cropped(:,:,sortingIndicies);

end