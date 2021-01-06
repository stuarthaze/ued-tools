function BeamBlockOut = fillBeamBlockEdges(EDGES,centre,nExpand)
% Input:
% EDGES - Logical 2d array
% centre - Structure containing centre.x and centre.y
% nExpand - integer, number of pixels to use to cover the edges
% Output:
% BeamBlockOut - logical 2d array of equal size to EDGES

cenLocation = round([centre.y,centre.x]); % input for imfill
FilledIMG = imfill(EDGES,cenLocation);
Block = double(FilledIMG - EDGES);
K = ones(nExpand);
BeamBlockExpanded = conv2(Block,K,'same');
BeamBlockOut = logical(BeamBlockExpanded);
end
