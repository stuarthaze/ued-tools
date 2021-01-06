function [IMGout,outliers] = outlierFilter2dMedian(IMG,threshold)
% Replaces pixels in IMGin with median value of neigbours if it differs by 
% more than threshold*std2(IMG-IMGmed)
    IMGmed = medfilt2(IMG);
    IMGdiff = IMG-IMGmed;
    stdIMG = std2(IMGdiff);
    outliers = (IMGdiff >= threshold*stdIMG) | (IMGdiff <= -threshold*stdIMG);
    IMGout = IMG;
    IMGout(outliers) = IMGmed(outliers);
end