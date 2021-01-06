function newCoords = moveMol(coords,sigma)
   nat = length(coords(:,1));
   dXYZ = sigma*randn(nat,3);
   newCoords = coords + dXYZ;
end