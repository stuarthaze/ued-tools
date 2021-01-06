function Iout = simpleScanProcessing(IMGstack,ROIstack)

nz = size(IMGstack,3);
nroi = size(ROIstack,3);
npix = sum(sum(ROIstack(:,:,1)))

Iout = zeros(nroi,nz);

for jj = 1:nz
    disp(['Processing img ',num2str(jj)]);
    for ii = 1:nroi
        Iout(ii,jj) = sum(sum(IMGstack(:,:,jj).*ROIstack(:,:,ii)))/npix;
    end
end
end