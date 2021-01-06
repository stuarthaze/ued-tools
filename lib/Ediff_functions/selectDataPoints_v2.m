function data2 = selectDataPoints_v2(data,tIndices,roiIndices,reliabilityCutoff)

data2 = data;
data2.spotNumbers = data.spotNumbers(roiIndices);
if isfield(data,'I_P') && ~isempty(data.I_P)
    data2.I_P = data.I_P(tIndices,roiIndices);
end
if isfield(data,'I_PP') && ~isempty(data.I_PP)
    data2.I_PP = data.I_PP(tIndices,roiIndices);
end

data2.ratios = data.ratios(tIndices,roiIndices);
data2.relativeChanges = data.relativeChanges(tIndices,roiIndices);
if isfield(data,'uncertainties_P') && ~isempty(data2.uncertainties_P)
    data2.uncertainties_P = data.uncertainties_P(tIndices,roiIndices);
end
if isfield(data,'uncertainties_PP') ~isempty(data2.uncertainties_P)
    data2.uncertainties_PP = data.uncertainties_PP(tIndices,roiIndices);
end
data2.uncertainties_R = data.uncertainties_R(tIndices,roiIndices);
data2.meanUncertainties_R = data.meanUncertainties_R(roiIndices);
data2.Wts = data.Wts(roiIndices,roiIndices);
disp('Warning - need to correct Wts for nConv > 1');
data2.HKL = data.HKL(roiIndices,:);
data2.nROIs = sum(roiIndices);
data2.timeDelays = data.timeDelays(tIndices);
data2.nTimeDelays = length(tIndices);
data2.initialResiduals = data.initialResiduals(tIndices,roiIndices);
data2.rmsResid_vs_T = data.rmsResid_vs_T(tIndices);
% Scale uncertainties with cutoff value 
% Note: This accounts for systematic error in intensities due to background
% or other factors based on discrepancy between theoretical and 
% experimental intensites of static pattern)
if nargin < 4
    reliabilityCutoff = 2;
end
relErrors = data.relativeErrors(roiIndices);
data2.relativeErrors = relErrors;
data2.meanUncertainties_R = data2.meanUncertainties_R./(1-(relErrors/reliabilityCutoff));
end