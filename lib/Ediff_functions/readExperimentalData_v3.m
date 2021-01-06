function FEDdata = readExperimentalData_v3(dirFED,dirIndexing,reliabilityCutoff)

addpath(dirFED);
addpath(dirIndexing);
load('Wts');
load('HKL');
load('Fobs_div_Fcalc');

reliabilityCutoffMin = 1/reliabilityCutoff;
spotIndices_reliable = (Fobs_div_Fcalc >= reliabilityCutoffMin) & (Fobs_div_Fcalc <= reliabilityCutoff);

% spotIndicies_reliable = spotVec((Fobs_div_Fcalc >= reliabilityCutoffMin) & (Fobs_div_Fcalc <= reliabilityCutoffMax));
% save('spotIndicies_reliable.mat','spotIndicies_reliable');

relativeChanges = load('relativeChanges.txt','-ascii');
FEDdata.relativeChanges = relativeChanges(:,spotIndices_reliable);
FEDdata.intensity_ratios = FEDdata.relativeChanges + 1;
uncertainties = load('relativeChanges_uncertainties.txt','-ascii');
FEDdata.uncertainties = uncertainties(:,spotIndices_reliable);
FEDdata.meanUncertainties = mean(FEDdata.uncertainties,1);  
FEDdata.Wts = Wts(:,spotIndices_reliable);
FEDdata.HKL = HKL(:,:,spotIndices_reliable);
nSpotsOrig = size(HKL,3);

Svec = 1:nSpotsOrig;
FEDdata.ROIs_originalNumbering = Svec(spotIndices_reliable);
FEDdata.nSpots = size(FEDdata.HKL,3);
FEDdata.nConv = size(FEDdata.HKL,1);
FEDdata.timeDelays = load('timeDelays.txt','-ascii');
FEDdata.Fobs_dif_Fcalc = Fobs_div_Fcalc;
FEDdata.spotIndices_reliable = spotIndices_reliable;
nTimeDelays = length(FEDdata.timeDelays);
FEDdata.nTimeDelays = nTimeDelays;
for T = 1:nTimeDelays
    initialResiduals(T,:) = FEDdata.relativeChanges(T,:)./FEDdata.meanUncertainties;
end
FEDdata.initialResiduals = initialResiduals;
FEDdata.rmsResid_vs_T = rms(initialResiduals,2);
end