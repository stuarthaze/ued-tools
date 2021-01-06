function IMGs = readSPE_images(DirectoryPath,string2match)
% Returns all images as a z-stack
files = dir(fullfile(DirectoryPath, string2match));
numFiles = length(files);
for indx = 1:numFiles
    IMGs(:,:,indx) = double(DataIOLibrary.DataIO.ReadSpe(fullfile(DirectoryPath,files(indx).name)));
end
end