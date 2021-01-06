function Mout = createDetectorCalibrationMat(sizeIMG,binning,filenameCalibration)

M1 = importdata(filenameCalibration);
sizeCalib = size(M1);

M2 = zeros(sizeIMG);

% calculate offset
halfsize = sizeIMG/2;
halfsizeCalib = sizeCalib/2;
offset = halfsizeCalib - halfsize*binning;

dX = offset(2) + (0:(sizeIMG(2)-1))*binning;
dY = offset(1) + (0:(sizeIMG(1)-1))*binning;

xBin = 1:binning;
yBin = 1:binning;

if binning == 1
    M2 = M1(dY+1,dX+1);
else
    for iX = 1:sizeIMG(2)
        for iY = 1:sizeIMG(1)
            X = dX(iX) + xBin;
            Y = dY(iY) + yBin;
            D = mean(mean(M1(Y,X)));
            M2(iY,iX) = D;
        end
    end
end
Mout = ones(size(M2))./M2;
end