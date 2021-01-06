function [IMG_2, num_good, IMG_stErrors] = fRemoveOutliersAndAverageIMG_IQR(IMGs,threshold)
% Takes a 3d array, or stack of images, IMGs(Y,X,img#)
% Removes pixels that differ from median by more than IQR*threshold
% IQR is the interquartile range which is close to 2 st.deviations for a
% Gaussian distr.
[ny,nx,nframes] = size(IMGs);
IMG_2 = zeros(ny,nx);
img_med = median(IMGs,3);
%img_mean = mean(IMGs,3);
img_m = img_med;
% img_std = std(IMGs,0,3);
IQR = quantile(IMGs,0.75,3)-quantile(IMGs,0.25,3);
img_uncert = IQR*threshold;
num_good = zeros(ny,nx);
for X = 1:nx
    for Y = 1:ny
        img_xy = IMGs(Y,X,:);
        good_pixels = (img_xy >= (img_m(Y,X)-img_uncert(Y,X))) ...
            & (img_xy <= (img_m(Y,X)+img_uncert(Y,X)));
        num_good(Y,X) = sum(good_pixels);
        values = img_xy(good_pixels);
        IMG_2(Y,X) = mean(values);
        IMG_stErrors(Y,X) = std(values)/sqrt(num_good(Y,X));
    end
end
% %check that std(IMGs) is not=0
% pixels_to_correct = (img_std == 0);
% IMG_2(pixels_to_correct) = img_mean(pixels_to_correct);
% num_good(pixels_to_correct) = nframes;
% IMG_stErrors(pixels_to_correct) = 0;