function IMG2 = viewTopNpeaks(M,locs,N,rCircle)
% SAH - Function to view peaks
% Input:    M - 2d array of peak heights
%           locs - 1d array of peak positions
%           N - number of peaks to view
%           rCircle - convlove peak with a circle (optional)
% Output:   IMG2 - 2d array with N peaks

IMG2 = zeros(size(M));
IMG2(locs(1:N)) = M(locs(1:N));
if nargin == 4
    IMG2 = convolveCircle(IMG2,rCircle);
end
end