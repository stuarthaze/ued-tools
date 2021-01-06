function RadialQuantile = fRadialQuantile(IMG,CEN,MASK,Q)
nx = size(IMG,2);
ny = size(IMG,1);
dmax(1) = sqrt(CEN.x^2+CEN.y^2);
dmax(2) = sqrt((CEN.x-nx)^2+CEN.y^2);
dmax(3) = sqrt(CEN.x^2+(ny-CEN.y)^2);
dmax(4) = sqrt((CEN.x-nx)^2+(ny-CEN.y)^2);
maxPix = round(max(dmax));
intensities{maxPix} = [];
for x = 1:nx
    for y = 1:ny
        if ~MASK(y,x)
            distSqrd = (x-CEN.x)^2 + (y-CEN.y)^2;
            r = round(sqrt(distSqrd));
            if distSqrd ~= 0;
                n = length(intensities{r});
                intensities{r}(n+1) = IMG(y,x);
            end
        end
    end
end
nr = length(intensities);
for R = 1:nr
    RadialQuantile(R) = quantile(intensities{R},Q);
end
