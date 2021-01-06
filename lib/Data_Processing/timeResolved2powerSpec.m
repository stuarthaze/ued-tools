function [PowerSpectrum, freqVec] = timeResolved2powerSpec(DATA,Tvec,nPad,plotResult)

nt = length(Tvec);
nROIs = size(DATA,2);
DATA2 = [DATA; zeros(nPad,nROIs)];

dataFFT = fft(DATA2);
dtVec = Tvec(2:nt)-Tvec(1:nt-1);
TvecMean = mean(Tvec);
if std(dtVec)/TvecMean > 1e-10
    disp('Warning: non-uniform step size');
end
dt1 = dtVec(1);
freqVec = calculateFrequenciesFFT(dt1,nt+nPad);
PowerSpectrum = dataFFT.*conj(dataFFT);

if plotResult
    nf_pos = ceil((nt+nPad)/2);
    figure(); 
    subplot(1,2,1); imagesc([],Tvec,DATA);
    title('TR data'); xlabel('ROI'); ylabel('Time / ps'); colormap('gray');
    subplot(1,2,2); imagesc([],freqVec(1:nf_pos),PowerSpectrum(1:nf_pos,:));
    title(['Power Spectrum, nData=',num2str(nt),', nPad=',num2str(nPad)]); xlabel('ROI'); ylabel('Frequency / THz');
    set(gca,'YDir','normal')
end

    function F = calculateFrequenciesFFT(dt,n)
    fs = 1/dt;
    df = fs/n;
    F0 = ((0:n-1)-floor(n/2))*df;
    F = ifftshift(F0);
    end

end