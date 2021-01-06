function radii = calcDistMat(nPixels,centre)
radii = zeros(nPixels.y,nPixels.x);
for X = 1:nPixels.x
    for Y = 1:nPixels.y
        radii(Y,X) = sqrt((X-centre.x)^2 + (Y-centre.y)^2);
    end
end
end