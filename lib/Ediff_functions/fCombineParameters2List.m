function P = fCombineParameters2List(X,U)
sizeX = size(X);
nX = sizeX(1)*sizeX(2);
sizeU = size(U);
nU = sizeU(1)*sizeU(2);
nParam = nX + nU;
P = zeros(nParam,1);
P(1:nX) = reshape(X,[nX,1]);
P((nX+1):nParam) = reshape(U,[nU,1]);
end

