function [IMGstack, timeDelays] = readSPE_TR(DirectoryPath,string2match)
dataIOAssembly = NET.addAssembly('\\win.desy.de\home\hayess\My Documents\MATLAB\lib\Data_Processing\DataIOLib.dll');
files = dir(fullfile(DirectoryPath, ['*',string2match,'*']));
numFiles = length(files);
for indx = 1:numFiles
    timeDelays(indx) = str2double(regexp(files(indx).name, '\d*(?=fs)', 'match'));
    IMGstack(:,:,indx) = double(DataIOLibrary.DataIO.ReadSpe(fullfile(DirectoryPath,files(indx).name)));
end
end