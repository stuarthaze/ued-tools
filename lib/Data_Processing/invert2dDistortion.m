function [X2,Y2] = invert2dDistortion(X1,Y1,nx,ny)
% S.Hayes Feb 2020
% Input: X1,Y1 - 2d matricies of x,y coordinates of distorted image
%        nx,ny - number of pixels for new coordinates
% Output: X2,Y2 - d2 matricies of coordinates in X1,Y1 corresponding to a
%                 linear range of coordinates 1:nx, 1:ny

X2 = ones([ny,nx]);
Y2 = ones([ny,nx]);

for x = 1:nx
    for y = 1:ny
        dx = X1 - x;
        dy = Y1 - y;
        r2 = dx.^2 + dy.^2;
        indx = find(r2 == min(r2,[],'all'));
        Y2(y,x) = mod(indx,ny);
        X2(y,x) = ceil(indx/ny);
    end
end