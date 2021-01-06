function [IMGs,timeDelays] = processScanAndor(dataDir,string2match,sortByTime)

files = dir(fullfile(dataDir, string2match));
num_imgs = length(files);
fname1 = files(1).name;
testimg = ReadspeAndor(fullfile(dataDir,fname1));
[NY,NX] = size(testimg); 
IMGs = zeros(NY,NX,num_imgs);
%IMG_BG = zeros(NY,NX);
for indx  = 1:num_imgs
    disp(['Reading file ',num2str(indx),' of ',num2str(num_imgs)]);
    fname = files(indx).name;
    IMGs(:,:,indx) = double(ReadspeAndor(fullfile(dataDir,fname)));
    
    searchPattern = '\d*(?=fs)';
    timeDelayCell = regexp(fname, searchPattern, 'match');
    timeDelays(indx) = str2double(timeDelayCell{1});
end

if (nargin == 3) && sortByTime
    [timeDelaysSorted, sortingIndicies] = sort(timeDelays);
    IMGsorted = IMGs(:,:,sortingIndicies);
    IMGs = IMGsorted;
    timeDelays = timeDelaysSorted;
end

end