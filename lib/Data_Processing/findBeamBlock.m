function [BeamBlockOut, EdgesExpanded] = findBeamBlock(IMG,centre,nExpand)
% Input:
%   IMG   - ordinary 2d matrix
%   centre - Structure containing centre.x and centre.y
%   nExpand - [int, int] number of pixels [y,x] to use to join the edges
% Output:
%   logical 2d array of equal size to IMG
if size(nExpand) == 1
    nExpand(2) = nExpand(1);
end
BW1 = edge(IMG,'sobel');
BW2 = edge(IMG,'canny');
EDGES = (BW1+BW2);
K1 = ones(nExpand);
EdgesExpanded = conv2(EDGES,K1,'same');
EdgesL = logical(EdgesExpanded); % imfill needs a logical
cenLocation = round([centre.y,centre.x]); % input for imfill
FilledIMG = imfill(EdgesL,cenLocation);
BeamBlock = FilledIMG - EdgesL;
nExpandBlock = nExpand+5;
K2 = ones(nExpandBlock);
% round off edges of convolution Kernal
K2(1,1) = 0;
K2(1,nExpandBlock(2)) = 0;
K2(nExpandBlock(1),1) = 0;
K2(nExpandBlock(1),nExpandBlock(2)) = 0;
BeamBlockExpanded = conv2(double(BeamBlock),K2,'same');
BeamBlockOut = logical(BeamBlockExpanded);

end