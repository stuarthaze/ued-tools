function data2 = selectDataPoints(data,tIndices,roiIndices)

data2 = data;
data2.spotNumbers = data.spotNumbers(roiIndices);
data2.I_P = data.I_P(tIndices,roiIndices);
data2.I_PP = data.I_PP(tIndices,roiIndices);

data2.ratios = data.ratios(tIndices,roiIndices);
data2.relativeChanges = data.relativeChanges(tIndices,roiIndices);
data2.uncertainties_P = data.uncertainties_P(tIndices,roiIndices);
data2.uncertainties_PP = data.uncertainties_PP(tIndices,roiIndices);
data2.uncertainties_R = data.uncertainties_R(tIndices,roiIndices);
data2.meanUncertainties_R = data.meanUncertainties_R(roiIndices);
data2.Wts = data.Wts(roiIndices,roiIndices);
disp('Warning - need to correct Wts for nConv > 1');
data2.HKL = data.HKL(roiIndices,:);
data2.nROIs = length(roiIndices);
data2.timeDelays = data.timeDelays(tIndices);
data2.nTimeDelays = length(tIndices);
data2.initialResiduals = data.initialResiduals(tIndices,roiIndices);
data2.rmsResid_vs_T = data.rmsResid_vs_T(tIndices);
end