function TR_data = processScan_medianIMGs(dataDir,ROIs3d,IMG_BG,imgCrop)

backupDirectories = dir(fullfile(dataDir,'Backup*'))
%Check number of frames per time point
directoryName = [dataDir, '\', backupDirectories(1).name];
filesPP = dir(fullfile(directoryName, '*Pump+Probe*'));
numFrames = length(filesPP);
NT = length(backupDirectories);
nPixSpot = sum(sum(ROIs3d(:,:,1)));
[NY_full,NX_full] = size(IMG_BG);
[NY_crop,NX_crop,Nspots] = size(ROIs3d);

IMGs_P = zeros(NY_full,NX_full,NT);
IMGs_PP = zeros(NY_full,NX_full,NT);

imgs_t_p = zeros(NY_full,NX_full,numFrames);
imgs_t_pp = zeros(NY_full,NX_full,numFrames);

spotIntensities_p = zeros(NT,Nspots,numFrames);
spotIntensities_pp = zeros(NT,Nspots,numFrames);

for tIndex = 1:NT
    fprintf(['Processing time point ',num2str(tIndex), ' of ', num2str(NT),'\n'])
    directoryName = [dataDir, '\', backupDirectories(tIndex).name];
    filesProbe = dir(fullfile(directoryName, '* Probe*'));
    filesPP = dir(fullfile(directoryName, '*Pump+Probe*'));
    searchPattern = '\d*(?=fs)';
    timeDelayCell = regexp(filesPP(1).name, searchPattern, 'match');
    timeDelays(tIndex) = str2double(timeDelayCell{1});
    for frame = 1:numFrames
        imgs_t_p(:,:,frame) = double(DataIOLibrary.DataIO.ReadSpe(fullfile(directoryName,filesProbe(frame).name)))-IMG_BG;
        imgs_t_pp(:,:,frame) = double(DataIOLibrary.DataIO.ReadSpe(fullfile(directoryName,filesPP(frame).name)))-IMG_BG;
        for spotIndex = 1:Nspots
            spotIntensities_p(tIndex,spotIndex,frame) = ...
                sum(sum(imgs_t_p(imgCrop.y,imgCrop.x,frame).*ROIs3d(:,:,spotIndex)))/nPixSpot;
            spotIntensities_pp(tIndex,spotIndex,frame) = ...
                sum(sum(imgs_t_pp(imgCrop.y,imgCrop.x,frame).*ROIs3d(:,:,spotIndex)))/nPixSpot;
        end
    end
    IMGs_P(:,:,tIndex) = median(imgs_t_p,3);
    IMGs_PP(:,:,tIndex) = median(imgs_t_pp,3);

end
[timeDelaysSorted, sortingIndicies] = sort(timeDelays);
IMG_diffs = IMGs_PP - IMGs_P;
TR_data.IMG_diffs = IMG_diffs(:,:,sortingIndicies);
TR_data.timeDelays = timeDelaysSorted;
TR_data.IMGs_P = IMGs_P(:,:,sortingIndicies);
TR_data.IMGs_PP = IMGs_PP(:,:,sortingIndicies);
TR_data.spotIntensities_P = spotIntensities_p(sortingIndicies,:,:);
TR_data.spotIntensities_PP = spotIntensities_pp(sortingIndicies,:,:);
end