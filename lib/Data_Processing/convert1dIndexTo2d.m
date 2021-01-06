function [m, n] = convert1dIndexTo2d(index,nM,nN)
n = ceil(index/nM);
m = index-(n-1)*nM;
end

