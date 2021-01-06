function img = xyz2image_withPartOccup(ATOMS,XYZ,dimension,xRange,yRange,AtmSize,OCC)
nX = length(xRange);
nY = length(yRange);
nAt = length(ATOMS);
img = zeros(nY,nX);
if dimension == 1
    a = 2;
    b = 3;
elseif dimension == 2
    a = 1;
    b = 3;
else
    a = 1;
    b = 2;
end
for x = 1:nX
    for y = 1:nY
        for atm = 1:nAt
            rAtSqrd = (XYZ(atm,a)-xRange(x))^2 + (XYZ(atm,b)-yRange(y))^2;
            img(y,x) = img(y,x) + OCC(atm)*ATOMS(atm)*exp(-rAtSqrd/(2*AtmSize^2));
        end
    end
end