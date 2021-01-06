function IMGs = readSPE_images_16bit(DirectoryPath,string2match)
% Returns all images as a z-stack
files = dir(fullfile(DirectoryPath, ['*',string2match,'*']));
numFiles = length(files);
fname = files(1).name;
testimg = double(DataIOLibrary.DataIO.ReadSpe(fullfile(DirectoryPath,fname)));
[NY,NX] = size(testimg); 
IMGs = zeros(NY,NX,numFiles,'uint16');

for indx = 1:numFiles
    IMGs(:,:,indx) = uint16(double(DataIOLibrary.DataIO.ReadSpe(fullfile(DirectoryPath,files(indx).name))));
end
end