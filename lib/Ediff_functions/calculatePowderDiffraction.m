function [IMG, S] = calculatePowderDiffraction(expParams,sampleParams)

ONES = ones(expParams.ny,expParams.nx);
[Xmat,Ymat] = meshgrid(1:expParams.nx, 1:expParams.ny);
xCenBeam = expParams.xCenBeam;
yCenBeam = expParams.yCenBeam;
sqrtPixelAspect = sqrt(expParams.pixelAspect);
[DX_lens, DY_lens] = recentreAndRescalePixels(Xmat,Ymat,expParams.xCenLens,expParams.yCenLens);
[X2_lens, Y2_lens] = radialLensDistortion(DX_lens,DY_lens);

[xCenBeamL1, yCenBeamL1] = recentreAndRescalePixels(xCenBeam,yCenBeam,expParams.xCenLens,expParams.yCenLens);
[xCenBeam2, yCenBeam2] = radialLensDistortion(xCenBeamL1,yCenBeamL1);


XX_beam = X2_lens - xCenBeam2;
YY_beam = Y2_lens - yCenBeam2;
[XX_beam2, YY_beam2] = astigmaticDistortion(XX_beam, YY_beam);
R = sqrt(XX_beam2.^2 + YY_beam2.^2);
S = 4*pi*sin(0.5*atan(R*expParams.pixelSize/expParams.distanceDet))/expParams.wavel;
BG = expParams.radialBackground(1)*exp(S/expParams.radialBackground(2));

ringWidthS = expParams.peakWidth*expParams.radialDist(1)*2*pi*expParams.pixelSize/(expParams.distanceDet*expParams.wavel);

IMG = expParams.intensityScale*simulateDiffractionRings(S,sampleParams.radii,sampleParams.intensities,ringWidthS) + BG;

    function [XX2,YY2] = radialLensDistortion(XX1,YY1)
    p = expParams.radialDist;
    R2rd = XX1.^2 + YY1.^2;
    RR = sqrt(R2rd);
    XX2 = XX1.*(p(1) + 0.001*p(2)*RR + 1e-6*p(3)*R2rd);
    YY2 = YY1.*(p(1) + 0.001*p(2)*RR + 1e-6*p(3)*R2rd);
    end

    function [XX2, YY2] = rotateCoordMatricies(XX1,YY1,theta)
    XX2 = XX1*cos(theta) - YY1*sin(theta);
    YY2 = XX1*sin(theta) + YY1*cos(theta);
    end

    function [XX2,YY2] = astigmaticDistortion(XX1,YY1)
    [XX2, YY2] = rotateCoordMatricies(XX1,YY1,expParams.thetaLensAstig);
    XX3 = XX2*expParams.lensAstigmatism;
    YY3 = YY2/expParams.lensAstigmatism;
    [XX2,YY2] = rotateCoordMatricies(XX3,YY3,-expParams.thetaLensAstig);
    end
    
    function [XX2, YY2] = recentreAndRescalePixels(XX1,YY1,xcen,ycen)
    XX2 = sqrtPixelAspect*(XX1 - xcen);
    YY2 = (YY1 - ycen)/sqrtPixelAspect;
    end

end