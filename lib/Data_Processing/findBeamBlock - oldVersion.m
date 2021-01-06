function BeamBlock = findBeamBlock(IMG,centre,xExpand)
% Input:
% IMG   - ordinary 2d matrix
% centre - Structure containing centre.x and centre.y
% xExpand - integer, number of pixels to use to cover the edges
% Output - logical 2d array of equal size to IMG
function IMG2 = expandEdges(IMG1)
    [ny,nx] = size(IMG1);
    Xshift(1,:) = (1:nx)+1;
    Xshift(1,nx) = nx;
    Xshift(2,:) = (1:nx)-1;
    Xshift(2,1) = 1;

    Yshift(1,:) = (1:ny)+1;
    Yshift(1,ny) = ny;
    Yshift(2,:) = (1:ny)-1;
    Yshift(2,1) = 1;

    IMG2 = IMG1 | IMG1(:,Xshift(1,:)) | IMG1(:,Xshift(2,:)) | IMG1(Yshift(1,:),:) | IMG1(Yshift(2,:),:);
end

BW1 = edge(IMG,'sobel');
BW2 = edge(IMG,'canny');
EDGES = (BW1+BW2);
EdgesL = logical(EDGES); % imfill needs a logical
Edges2 = expandEdges(EdgesL);
cenLocation = round([centre.y,centre.x]); % input for imfill
FilledIMG = imfill(Edges2,cenLocation);
BeamBlock = FilledIMG - Edges2;
for N = 1:xExpand
    BeamBlock = expandEdges(BeamBlock);
end

end