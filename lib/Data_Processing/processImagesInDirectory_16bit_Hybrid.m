function [IMG_out, IMG_number_good_hybrid]  = processImagesInDirectory_16bit_Hybrid(dataDir, string2match, rejectionThreshold, ratioNormQuant, displayIMG, displayNumberRejected)
% Reads images in folder containing [string2match],
% removes outliers from z-stack and averages remaining pixels
% Reads images as unit16 datatypes - only use for raw images 
% (not sums or differences)

% dataLocation;
files = dir(fullfile(dataDir, string2match));
num_imgs = length(files);
fname = files(1).name;
testimg = double(DataIOLibrary.DataIO.ReadSpe(fullfile(dataDir,fname)));
[NY,NX] = size(testimg); 
IMGs = zeros(NY,NX,num_imgs,'uint16');
%IMG_BG = zeros(NY,NX);
disp('Reading Data');
for indx  = 1:num_imgs
    fname = files(indx).name;
    IMGs(:,:,indx) = uint16(double(DataIOLibrary.DataIO.ReadSpe(fullfile(dataDir,fname))));
end

disp('Removing Outliers');
[IMG_out, IMG_number_good_hybrid] = fRemoveOutliersAndAverageIMG_16bitIMGstack_Hybrid(IMGs,rejectionThreshold,ratioNormQuant);
% Estimates value as:                   R*mean + (1-R)*median
% Estimates uncertainty as:             R*std + (1-R)*IQR*iqr2std
% Removes pixels outside the range:     estimate +- uncertainty*threshold
if displayIMG
    figure(); imagesc(IMG_out);
    colorbar();
    title('Averaged IMG');
end

if displayNumberRejected
    figure(); imagesc(num_imgs-IMG_number_good_hybrid);
    colorbar();
    title('number of background rejections (hybrid alg.)');
end
end