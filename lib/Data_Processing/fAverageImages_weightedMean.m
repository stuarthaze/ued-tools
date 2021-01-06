function [IMG_2, IMG_var] = fAverageImages_weightedMean(IMGs)
% Takes a 3d array, or stack of images, IMGs(Y,X,img#)
% Calculates mean & variance for each pixel
% Calculates weighted mean: wt = exp(dx^2/var)
[ny,nx,nz] = size(IMGs);
img_mean = mean(IMGs,3);
img_var = var(IMGs,0,3);
dI_3d = zeros([ny,nx,nz]);
wts_3d = zeros([ny,nx,nz]);

for Z = 1:nz
    dI_2d = IMGs(:,:,Z) - img_mean;
    wts_3d(:,:,Z) = exp(-0.5*dI_2d.^2./img_var);
end

wts_sum = sum(wts_3d,3);
IMGs_wtd = IMGs.*wts_3d;
IMG_2 = sum(IMGs_wtd,3)./wts_sum;

for Z = 1:nz
    dI_3d(:,:,Z) = IMGs(:,:,Z) - IMG_2;
end

IMG_var = sum(wts_3d.^2 * dI_3d.^2,3)./(wts_sum.^2);

end
