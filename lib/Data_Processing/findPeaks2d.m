function PKS = findPeaks2d(IMG,Width)
% SAH 20.06.2019
% Input:    IMG - 2d array
%           Width - Window size for peak search in pixels
% Output:   PKS - 2d array with non-zero values of peak heights at peak
%                   positions


[ny,nx] = size(IMG);
PKS = zeros(ny,nx);
% Reduce width to nearest odd int for a symmetric window
HW = floor((Width-1)/2); 
W = HW*2;

for Y = 1:(ny-W)
    for X = 1:(nx-W)
        Ipixel = IMG(Y+HW,X+HW);
        I_window = IMG(Y:(Y+W), X:(X+W));
        if Ipixel == max(max(I_window))
            PKS(Y+HW,X+HW) = Ipixel;
        end
    end
end

end
