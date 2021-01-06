function [IMG_2, num_good, IMG_stErrors] = fRemoveOutliersAndAverageIMG_IQR_16bitIMGstack(IMGs,threshold)
% Takes a 3d array, or stack of images, IMGs(Y,X,img#)
% Removes pixels that differ from median by more than IQR*threshold
% IQR is the interquartile range which is close to 2 st.deviations for a
% Gaussian distr.
[ny,nx,nframes] = size(IMGs);
IMG_2 = zeros(ny,nx);
num_good = zeros(ny,nx);
IMG_stErrors = zeros(ny,nx);
Z = zeros(1,nframes);
good_pixels = logical(Z);
for Y = 1:ny
    disp(['Processing row ', num2str(Y),'\n'])
    for X = 1:nx
        Z = double(IMGs(Y,X,:));
        med_z = median(Z);
        iqr_z = iqr(Z);
        uncert = threshold*iqr_z;
        good_pixels = (Z >= (med_z-uncert)) & (Z <= (med_z+uncert));
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