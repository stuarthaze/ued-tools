function dx = calcNearestNeighbourDiffs(x)

nx = length(x);
dx = x(2:nx) - x(1:nx-1);
end