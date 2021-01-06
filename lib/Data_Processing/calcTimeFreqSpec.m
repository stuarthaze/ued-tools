function [FTcomplex, frequencies] = calcTimeFreqSpec(DATA,dt)

nData = length(DATA);
nFreqs = floor(nData/2);
FTcomplex = complex(zeros(nFreqs,nData));
xvec = 1:nData;
freqs = (1:nFreqs)*0.5/nFreqs;
windowSizes = 1./freqs;
frequencies = freqs/dt;

for t = 1:nData
    for findx = 1:nFreqs
        window = exp(-(xvec-t).^2/windowSizes(findx)^2)/windowSizes(findx);
        Wavelet = window.*exp(2i*pi*freqs(findx)*(xvec-t));
        FTcomplex(findx,t) = sum(Wavelet.*DATA');
    end
end
%     
% 
%     function GW = gaussianWindow(x0,sigma)
%         GW = exp(-0.5*(xvec-x0).^2/sigma^2)/sigma;
%     end

end