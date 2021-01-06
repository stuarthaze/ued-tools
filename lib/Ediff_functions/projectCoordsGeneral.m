function IMG = projectCoordsGeneral(XYZ,atoms,ViewVector,VerticalDirection,X,Y,atomSize)
% XYZ    - Cartesian Coordinates
% ViewVector - Direction of view
% VerticalDirection - Vector used to determine the vertical direction
% X, Y - vectors of coordinates to use in the image, IMGxIm
% atomSize - size of atoms in IMG
% IMG - 2D projection of XYZ coordinates
%--------------------------------------------
% First define new coordinates from the point of view - ViewVector = Z*a
X2 = cross(ViewVector,VerticalDirection);
Y2 = -cross(ViewVector,X2);
XYZ2_basis(1,:) = X2/sqrt(sum(X2.^2));
XYZ2_basis(2,:) = Y2/sqrt(sum(Y2.^2));
XYZ2_basis(3,:) = ViewVector/sqrt(ViewVector.^2);
%-------------------------------------
% Transform coordinates
nAt = size(XYZ,1);
XYZ2 = XYZ*(XYZ2_basis');
%------------------------------------------
% Create IMG
nX = length(X);
nY = length(Y);
IMG = zeros(nY,nX);
for a = 1:nY
    for b = 1:nX
        for atm = 1:nAt
            dX = X(b) - XYZ2(atm,1);
            dY = Y(1-a+nY) - XYZ2(atm,2);
            IMG(a,b) = IMG(a,b)+atoms(atm)*exp(-(dX^2 + dY^2)/(atomSize^2));
        end
    end
end
end
