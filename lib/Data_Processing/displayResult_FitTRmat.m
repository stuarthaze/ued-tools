function displayResult_FitTRmat(fitStruct)

data = fitStruct.Mdata;
tDelays_zeroCorr = fitStruct.timeDelays - fitStruct.t0;
ndat = size(data,1)*size(data,2);
p=0.002;
scaleMax = quantile(reshape(data,[ndat,1]),1-p);
scaleMin = quantile(reshape(data,[ndat,1]),p);
SCALE = [scaleMin,scaleMax];
figure(); 
subplot(1,3,1); imagesc(data,SCALE); title('Data');
subplot(1,3,2); imagesc(fitStruct.Mcalc,SCALE); title(['Fitted, \sigma=', num2str(round(fitStruct.sigmaIRF)),' fs']);
subplot(1,3,3); imagesc(fitStruct.Mdiff,SCALE ); colormap('gray'); title('Residuals');
end