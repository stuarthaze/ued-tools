function ROIstack = positions2roiMasks(X,Y,size,diam)

nROIs = length(X);
radius = diam*0.5;
ROIstack = zeros([size.y,size.x,nROIs]);
[XX,YY] = meshgrid(1:size.x, 1:size.y);

for ii = 1:nROIs
    dXX = XX - X(ii);
    dYY = YY - Y(ii);
    dRR = sqrt(dXX.^2 + dYY.^2);
    ROIstack(:,:,ii) = dRR <= radius;
end
end