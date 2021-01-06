function rAB = atomPairDistr(atA,atB,ensemble)
   nMols = length(ensemble(1,1,:));
   nAts = length(ensemble(:,1,1));
   rAB = zeros(1,nMols);
   rAB2 = zeros(1,3);
   for m = 1:nMols
       rAB2(:) = (ensemble(atA,:,m)-ensemble(atB,:,m)).^2;
       rAB(m) = sqrt(sum(rAB2));
   end
end