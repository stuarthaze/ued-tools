function IMG2 = maskCentre(IMG,centre,radius)

[ny,nx] = size(IMG);
IMG2 = zeros(size(IMG));
[XX,YY] = meshgrid(1:nx,1:ny);
YY = YY - centre.y;
XX = XX - centre.x;
RR2rd = YY.^2 + XX.^2;
MASK = RR2rd > (radius^2);
IMG2 = IMG.*MASK;
end