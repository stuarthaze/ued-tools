function [BG_2d, Radial_BG] = radialQuantile_BG(IMG, centre, MASK, quantile)

Radial_quantile = fRadialQuantile(IMG,centre,MASK,quantile);
R = 1:length(Radial_quantile);
RRadial = R.*Radial_quantile;
RRadial_NaNedit = replaceNaNsWithNearestValue(RRadial);
RRadial_smoothed = fGaussianFilterVertical(RRadial_NaNedit',5)';
Radial_BG = RRadial_smoothed./R;
IMG_Rmat = IMG_calcDists2Centre(IMG,centre);
IMG_Rmat(IMG_Rmat < 1) = 1;
BG_2d = Radial_BG(round(IMG_Rmat));

