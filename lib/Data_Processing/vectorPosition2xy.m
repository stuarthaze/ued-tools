function [x,y] = vectorPosition2xy(position,nx,ny)
% returns the [x,y] coordinates from a position in a vector
% upon transforming to a 2d array with nx columns and ny rows
x = ceil(position/ny);
y = mod(position,ny);
y(find(y==0)) = ny;
end