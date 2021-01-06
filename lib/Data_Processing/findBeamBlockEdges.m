function EdgesLogical = findBeamBlockEdges(IMG, nExpand, showResult)
% Input:
% IMG   - ordinary 2d matrix
% centre - Structure containing centre.x and centre.y
% nExpand - integer, size of convolution kernal for expansion
% Output - logical 2d array of equal size to IMG
%--------------------------------------------------
BW1 = edge(IMG,'sobel');
BW2 = edge(IMG,'canny');
EDGES = BW1+BW2;
Kernal = ones(nExpand);
EdgesExpanded = conv2(EDGES,Kernal,'same');
EdgesLogical = logical(EdgesExpanded);
if showResult
    figure(); imagesc(Edges2); pbaspect([1 1 1]);
end
end