function par = expStruc2params(experimentStruct)

par(1) = experimentStruct.xCenBeam;
par(2) = experimentStruct.yCenBeam;
par(3) = experimentStruct.xCenLens;
par(4) = experimentStruct.yCenLens;
par(5:7) = experimentStruct.radialDist;
par(8) = experimentStruct.intensityScale;
par(9) = experimentStruct.lensAstigmatism;
par(10) = experimentStruct.thetaLensAstig;
par(11) = experimentStruct.peakWidth;
par([12, 13]) = experimentStruct.radialBackground;

end