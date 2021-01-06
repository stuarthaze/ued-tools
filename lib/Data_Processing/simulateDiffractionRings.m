function IMG = simulateDiffractionRings(Rmat,radii,intensities,ringWidth)

nRings = length(radii);
[ny,nx] = size(Rmat);
IMGstack = zeros([ny,nx,nRings]);

for indx = 1:nRings
    dR = Rmat-radii(indx);
    IMGstack(:,:,indx) = intensities(indx)*exp(-0.5*dR.^2/(ringWidth^2));
end
IMG = sum(IMGstack,3);
end
