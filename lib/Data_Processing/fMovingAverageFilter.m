function IMG2 = fMovingAverageFilter(IMG,nPixels)

nPix = round(nPixels);
K = ones(nPix)/nPix^2;  
IMG2 = conv2(IMG,K,'same');
end