function data2 = selectTimePoints(data,tIndicies)

data2 = data;
data2.relativeChanges = data.relativeChanges(tIndicies,:);
data2.ratios = data.ratios(tIndicies,:);
data2.uncertainties_R = data.uncertainties_R(tIndicies,:);
data2.timeDelays = data.timeDelays(tIndicies);
data2.nTimeDelays = length(tIndicies);
data2.initialResiduals = data.initialResiduals(tIndicies,:);
end