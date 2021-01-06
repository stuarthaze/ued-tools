function img = readSPE( dataDir,filename )
% dataIOAssembly = NET.addAssembly('\\win.desy.de\home\hayess\My Documents\MATLAB\lib\Data_Processing\DataIOLib.dll');
img = DataIOLibrary.DataIO.ReadSpe(fullfile(dataDir,filename)); 
img = double(img);
end
