function IMGs2 = binImages3d(IMGs,binning)
% Bins a z-stack of images

[ny1,nx1,nz] = size(IMGs);
nx2 = floor(nx1/binning);
ny2 = floor(ny1/binning);
dx = 1:binning;
IMGs2 = zeros([ny2,nx2,nz]);
for x2 = 1:nx2
    for y2 = 1:ny2
        x1 = binning*(x2-1) + dx;
        y1 = binning*(y2-1) + dx;
        IMGs2(y2,x2,:) = mean(mean(IMGs(y1,x1,:),2),1);
    end
end

end