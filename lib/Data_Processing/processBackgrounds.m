function IMG_bg = processBackgrounds(dataDir, rejectionThreshold, displayIMG, displayNumberRejected)
% Reads contents of folder containing backgrounds removes outliers and 
% averages

% dataLocation;
files = dir(fullfile(dataDir, '*img*'));
num_bgs = length(files);
fname = files(1).name;
testimg = double(DataIOLibrary.DataIO.ReadSpe(fullfile(dataDir,fname)));
[NY,NX] = size(testimg); 
BGs = zeros(NY,NX,num_bgs);
%IMG_BG = zeros(NY,NX);
for indx  = 1:num_bgs
    fname = files(indx).name;
    BGs(:,:,indx) = double(DataIOLibrary.DataIO.ReadSpe(fullfile(dataDir,fname)));
end

[IMG_bg, IMG_bg_number_good_hybrid] = fRemoveOutliersAndAverageIMG_hybrid(BGs,rejectionThreshold);
if displayIMG
    figure(); imagesc(IMG_bg);
    colorbar();
    title('Averaged IMG');
end

if displayNumberRejected
    figure(); imagesc(num_bgs-IMG_bg_number_good_hybrid);
    colorbar();
    title('number of background rejections (hybrid alg.)');
end
end