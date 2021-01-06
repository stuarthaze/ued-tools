function [XYZ2, Tmat] = transformCoords_fromViewVector(XYZ,ViewVector,VerticalDirection)
% XYZ    - Cartesian Coordinates
% ViewVector - Direction of view (new Z)
% VerticalDirection - Vector used to determine the vertical direction (coplanar with new Y)
%--------------------------------------------
% First define new coordinates from the point of view - ViewVector = Z*a
X2 = -cross(ViewVector,VerticalDirection);
Y2 = cross(ViewVector,X2);
XYZ2_basis(1,:) = X2/norm(X2);
XYZ2_basis(2,:) = Y2/norm(Y2);
XYZ2_basis(3,:) = ViewVector/norm(ViewVector);
Tmat = inv(XYZ2_basis);
%-------------------------------------
% Transform coordinates
XYZ2 = XYZ*Tmat;