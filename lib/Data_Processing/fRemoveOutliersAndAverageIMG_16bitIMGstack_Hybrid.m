function [IMG_2, num_good, IMG_stErrors] = fRemoveOutliersAndAverageIMG_16bitIMGstack_Hybrid(IMGs,threshold,ratioNormalQuantile)
% Takes a 3d array, or stack of images, IMGs(Y,X,img#)
% Estimates value as:                   R*mean + (1-R)*median
% Estimates uncertainty as:             R*std + (1-R)*IQR*iqr2std
% Removes pixels outside the range:     estimate +- uncertainty*threshold

iqr2sigma = 1.349;
[ny,nx,nframes] = size(IMGs);
IMG_2 = zeros(ny,nx);
num_good = zeros(ny,nx);
IMG_stErrors = zeros(ny,nx);
Z = zeros(1,nframes);
good_pixels = logical(Z);
for Y = 1:ny
%     disp(['Processing row ', num2str(Y)])
    for X = 1:nx
        Z = double(IMGs(Y,X,:));
        estimate_z = ratioNormalQuantile*mean(Z) + (1-ratioNormalQuantile)*median(Z);
        uncertainty = ratioNormalQuantile*std(Z) + (1-ratioNormalQuantile)*iqr(Z)*iqr2sigma;
        good_pixels = (Z >= (estimate_z-uncertainty*threshold)) & (Z <= (estimate_z+uncertainty*threshold));
        num_good(Y,X) = sum(good_pixels);
        values = Z(good_pixels);
        IMG_2(Y,X) = mean(values);
        IMG_stErrors(Y,X) = std(values)/sqrt(num_good(Y,X));
    end
end
end
% %check that std(IMGs) is not=0
% pixels_to_correct = (img_std == 0);
% IMG_2(pixels_to_correct) = img_mean(pixels_to_correct);
% num_good(pixels_to_correct) = nframes;
% IMG_stErrors(pixels_to_correct) = 0;