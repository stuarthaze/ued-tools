function PKS = findPeaks2d_maxInWindow(IMG,Width,threshold)

[ny,nx] = size(IMG);
PKS = zeros(ny,nx);
W= Width-1;
L = Width^2;
HW = floor(Width/2);

for Y = 1:(ny-W)
    Y
    for X = 1:(nx-W)
        if IMG(Y+HW,X+HW) >= threshold
            V = reshape(IMG(Y:(Y+W), X:(X+W)),1,L);
            [val,indx] = max(V);
            YX = convert1dIndexTo2d(indx,Width,Width);
            YX2 = YX + [Y,X] - [1,1];
            PKS(YX2(1),YX2(2)) = PKS(YX2(1),YX2(2)) + val;
        end
    end
end

function P = convert1dIndexTo2d(index,nM,nN)
P = zeros(1,2);
n = ceil(index/nM);
m = index-(n-1)*nM;
P = [m,n];
end

end
