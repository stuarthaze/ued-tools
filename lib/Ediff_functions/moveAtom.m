function newCoords = moveAtom(Atm,coords,sigma)
   dXYZ = random('norm',0,sigma,[1,3]);
   newCoords = coords;
   newCoords(Atm,:) = coords(Atm,:) + dXYZ;
end