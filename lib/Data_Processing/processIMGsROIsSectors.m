function DATAout = processIMGsROIsSectors(IMGstack,ROIs,ROIsectors)


nROIs = size(ROIs,3);
nIMGs = size(IMGstack,3);

if nargin == 2

    npixROI = reshape(sum(sum(ROIs,2),1),[1,nROIs]);
    DATAout = zeros([nIMGs,nROIs]);

    for img = 1:nIMGs
        for roi = 1:nROIs
            DATAout(img,roi) = sum(sum(ROIs(:,:,roi).*IMGstack(:,:,img)))/npixROI(roi);
        end
    end
    
elseif nargin == 3
    nSect = size(ROIsectors,3);
    
    DATAout = zeros([nIMGs,nROIs,nSect]);

    for img = 1:nIMGs
        IMG = IMGstack(:,:,img);
        for roi = 1:nROIs
            for sec = 1:nSect
                ROI = ROIs(:,:,roi) & ROIsectors(:,:,sec);
                DATAout(img,roi,sec) = mean(IMG(ROI));
            end
        end
    end
else
    disp('error in number of input parameters');
end
    
end