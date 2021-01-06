function xyz2 = generateCoordsByInversion(xyz1,inversionCentre)
% Generates new coordinates, xyz2, by applying an inversion 
% operation to xyz1
nAt = size(xyz1,1);
moveXYZ = zeros(nAt,3);
for dim = 1:3
    moveXYZ(:,dim) = 2*inversionCentre(dim);
end
xyz2 = -xyz1 + moveXYZ;
end