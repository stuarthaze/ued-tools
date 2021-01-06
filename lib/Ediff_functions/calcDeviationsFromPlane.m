function dXplane = calcDeviationsFromPlane(XYZ,AtomsInPlane)
xyz = XYZ(AtomsInPlane,:);
nAt = size(xyz,1);
[R0 N] = calcPlane(xyz);
dXplane = zeros(1,nAt);
for atm = 1:nAt
    dxyz = xyz(atm,:)-R0;
    dXplane(atm) = dot(dxyz,N);
end
    function [r0 n] = calcPlane(X)
       r0 = mean(X);
       nAtX = size(X,1);
       nPairs = 0;
       for a = 1:(nAtX-1)
           for b = (a+1):nAtX
               nPairs = nPairs+1;
               Ra = X(a,:)-r0;
               Rb = X(b,:)-r0;
               vN(nPairs,:) = cross(Ra,Rb);
           end
       end
       refVect = vN(1,:);
       for vect = 1:nPairs
           direction = dot(vN(vect,:),refVect);
           if direction < 0
               vN(vect,:) = -vN(vect,:);
           end
       end
       V = sum(vN);
       n = V/norm(V);
   end
end
