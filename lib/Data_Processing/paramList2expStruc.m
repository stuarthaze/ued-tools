function experimentStruct = paramList2expStruc(par,expStruct0)

experimentStruct = expStruct0;
experimentStruct.xCenBeam = par(1);
experimentStruct.yCenBeam = par(2);
experimentStruct.xCenLens = par(3);
experimentStruct.yCenLens = par(4);
experimentStruct.radialDist = par(5:7);
experimentStruct.intensityScale = par(8);
experimentStruct.lensAstigmatism = par(9);
experimentStruct.thetaLensAstig = par(10);
experimentStruct.peakWidth = par(11);
experimentStruct.radialBackground = par([12, 13]);

end