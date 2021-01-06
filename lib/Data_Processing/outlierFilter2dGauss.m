function [IMGout,outliers] = outlierFilter2dGauss(IMG,threshold)
% Replaces pixels in IMGin with median value of neigbours if it differs by 
% more than threshold*std2(IMG-IMGmed)
    IMGfilt = imgaussfilt(IMG,0.5);
    IMGdiff = IMG-IMGfilt;
    stdIMG = std2(IMGdiff);
    outliers = (IMGdiff >= threshold*stdIMG) | (IMGdiff <= -threshold*stdIMG);
    IMGout = IMG;
    IMGout(outliers) = IMGfilt(outliers);
end