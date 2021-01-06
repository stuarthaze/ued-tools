function DATAout = processIMGsROIs(IMGstack,ROIs)

nROIs = size(ROIs,3);
nIMGs = size(IMGstack,3);
npixROI = reshape(sum(sum(ROIs,2),1),[1,nROIs]);
DATAout = zeros([nIMGs,nROIs]);

for img = 1:nIMGs
    for roi = 1:nROIs
        DATAout(img,roi) = sum(sum(ROIs(:,:,roi).*IMGstack(:,:,img)))/npixROI(roi);
    end
end
