function angles = calcBondAngles(XYZ,atomsABC)
xyzAt1 = XYZ(atomsABC(:,1),:);
xyzAt2 = XYZ(atomsABC(:,2),:);
xyzAt3 = XYZ(atomsABC(:,3),:);
bondVec1 = xyzAt1 - xyzAt2;
bondVec2 = xyzAt3 - xyzAt2;
% bondLengths1 = sum(bondVec1.^2,2)
normalization = sqrt(sum(bondVec1.^2,2)).*sqrt(sum(bondVec2.^2,2));
dotProd = dot(bondVec1,bondVec2,2)./normalization;
angles = acos(dotProd)';
end