function IMGsOut = outlierFilterGeneral(IMGs,filterOptions)
% IMGs: simgle image or stack of images
% options.type = 'Gaussian' or 'Median'
% options.window: Size of window (Median) or number of sigmas.
% options.ndim = 2 or 3 - number of dimensions of filter
% options.thresholdType = 'rms' or 'absolute'
% options.threshold = threshold value - multiplied by approximate noise level if thresholdType = 'rms'

    ndimIMG = length(size(IMGs));
    ndimFilt = filterOptions.ndim;
    if ndimIMG <= ndimFilt
        IMGsOut = filterIMG(IMGs);
    else
        imgSize = size(IMGs);
        IMGsOut = zeros(imgSize);
        for z = 1:imgSize(3)
            IMGsOut(:,:,z) = filterIMG(IMGs(:,:,z));
        end
    end

    function imgB = filterIMG(imgA)
        if strcmp(filterOptions.type,'Gaussian')
            if ndimFilt == 2
                imgF = imgaussfilt(imgA,filterOptions.window);
            elseif ndimFilt == 3
                imgF = imgaussfilt3(imgA,filterOptions.window);
            end
        elseif strcmp(filterOptions.type,'Median')
            if ndimFilt == 2
                window = filterOptions.window;
                if length(filterOptions.window) == 1
                    window = [filterOptions.window,filterOptions.window];
                end
                imgF = medfilt2(imgA,window);
            elseif ndimFilt == 3
                imgF = medfilt3(imgA,filterOptions.window);
            end
        else
            disp('Specify filter type');
            return
        end
        imgDiff = imgA - imgF;
        if strcmp(filterOptions.thresholdType,'rms')
            [nya,nxa,nza] = size(imgA);
            npix = nxa*nya*nza;
            threshold = rms(reshape(imgDiff,[1,npix]))*filterOptions.threshold;
        elseif strcmp(filterOptions.thresholdType,'absolute')
            threshold = filterOptions.threshold;
        end
            
        outliers = (imgDiff >= threshold) | (imgDiff <= -threshold);
        imgB = imgA;
        imgB(outliers) = imgF(outliers);
    end
end