function newCoords = moveMol(coords,sigma)
   nat = length(coords(:,1));
   dXYZ = random('norm',0,sigma,[nat,3]);
   newCoords = coords + dXYZ;
end