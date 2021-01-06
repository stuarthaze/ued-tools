function A = createNearestNeighbourDiffMatrix(dimensions)
% Creates matrix A with input dimensions [y, x]
A = zeros(dimensions);
iEnd = min(dimensions);
for ii = 1:iEnd
    A(ii,ii) = 1;
end
for ii = 1:(iEnd-1)
    A(ii+1,ii) = -0.5;
    A(ii,ii+1) = -0.5;
end
A(1,1) = 0.5;
A(iEnd,iEnd) = 0.5;

end

