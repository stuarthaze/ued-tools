function Itot = radialAverage(IMG,MASK,xCen,yCen)
imgSize = size(IMG);
nx = imgSize(2);
ny = imgSize(1);
rMax(1) = round(sqrt((xCen - 1)^2 + (yCen-1)^2));
rMax(2) = round(sqrt((xCen - nx)^2 + (yCen-1)^2));
rMax(3) = round(sqrt((xCen - 1)^2 + (yCen-ny)^2));
rMax(4) = round(sqrt((xCen - nx)^2 + (yCen-ny)^2));
maxDist = max(rMax);
Isum = zeros(2,maxDist);

for y = 1:ny
    for x = 1:nx
        if MASK(y,x)
            R = sqrt((y-yCen)^2 + (x-xCen)^2);
            nR = round(R);
            Isum(1,nR) = Isum(1,nR) + IMG(y,x);
            Isum(2,nR) = Isum(2,nR) + 1;
        end
    end
end

Itot = Isum(1,:)./Isum(2,:);
