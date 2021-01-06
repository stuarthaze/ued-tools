function IMG2 = windowFilter( IMG, minVal, maxVal )
%Filters an image, retaining all values greater than minVal and less
% than maxVal, setting all others to zero, producing a new image, IMG2
mask = (IMG > minVal) & (IMG < maxVal);
IMG2 = IMG.*mask;
end