function IMG2 = expandEdgesIter(IMG1,nIterations)

IMGL = logical(IMG1);
function img2 = expandEdges(img1)
    [ny,nx] = size(img1);
    Xshift(1,:) = (1:nx)+1;
    Xshift(1,nx) = nx;
    Xshift(2,:) = (1:nx)-1;
    Xshift(2,1) = 1;

    Yshift(1,:) = (1:ny)+1;
    Yshift(1,ny) = ny;
    Yshift(2,:) = (1:ny)-1;
    Yshift(2,1) = 1;

    img2 = img1 | img1(:,Xshift(1,:)) | img1(:,Xshift(2,:)) | img1(Yshift(1,:),:) | img1(Yshift(2,:),:);
end

for N = 1:nIterations
    IMGL = expandEdges(IMGL);
end

IMG2 = IMGL;
end

