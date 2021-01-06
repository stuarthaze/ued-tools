function [FTcomplex, FTpower] = calcSTFTcosWindow(DATA)

nData = length(DATA);
FTcomplex = complex(zeros(nData));
xvec = 1:nData;
fvec = 1./(1:nData);
sigmaVec = 1:nData;

for t = 1:nData
    for findx = 1:nData
        W = exp(-0.5*(xvec-t).^2/sigmaVec(findx)^2)/sigmaVec(findx).*exp(2i*pi*fvec(findx)*(xvec-t));
        FTcomplex(t,findx) = sum(W.*DATA');
    end
end
%     
% 
%     function GW = gaussianWindow(x0,sigma)
%         GW = exp(-0.5*(xvec-x0).^2/sigma^2)/sigma;
%     end

end