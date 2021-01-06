function R_IMG = IMG_calcDists2Centre(testIMG,centre)
[ny,nx] = size(testIMG);
Y = (1:ny)-centre.y;
X = (1:nx)-centre.x;
[XX,YY] = meshgrid(X,Y);
R_IMG = sqrt(XX.^2 + YY.^2);
end