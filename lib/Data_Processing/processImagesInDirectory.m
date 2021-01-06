function IMG_out = processImagesInDirectory(dataDir, string2match, rejectionThreshold, displayIMG, displayNumberRejected)
% Reads images in folder containing [string2match],
% removes outliers from z-stack and averages remaining pixels

% dataLocation;
files = dir(fullfile(dataDir, string2match));
num_imgs = length(files);
fname = files(1).name;
testimg = double(DataIOLibrary.DataIO.ReadSpe(fullfile(dataDir,fname)));
[NY,NX] = size(testimg); 
IMGs = zeros(NY,NX,num_imgs);
%IMG_BG = zeros(NY,NX);
for indx  = 1:num_imgs
    fname = files(indx).name;
    IMGs(:,:,indx) = double(DataIOLibrary.DataIO.ReadSpe(fullfile(dataDir,fname)));
end

[IMG_out, IMG_number_good_hybrid] = fRemoveOutliersAndAverageIMG_IQR(IMGs,rejectionThreshold);
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