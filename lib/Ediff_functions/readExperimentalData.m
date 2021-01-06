function FEDdata = readExperimentalData(dirFED,dirIndexing)

addpath(dirFED);
addpath(dirIndexing);
load('Wts');
load('HKL');
load('spotIndicies_reliable');
relativeChanges = load('relativeChanges.txt','-ascii');
FEDdata.relativeChanges = relativeChanges(:,spotIndicies_reliable);
FEDdata.intensity_ratios = FEDdata.relativeChanges + 1;
uncertainties = load('relativeChanges_uncertainties.txt','-ascii');
FEDdata.uncertainties = uncertainties(:,spotIndicies_reliable);
FEDdata.meanUncertainties = mean(FEDdata.uncertainties,1);  
FEDdata.Wts = Wts(:,spotIndicies_reliable);
FEDdata.HKL = HKL(:,:,spotIndicies_reliable);
FEDdata.nSpots = size(FEDdata.HKL,3);
FEDdata.nConv = size(FEDdata.HKL,1);
FEDdata.timeDelays = load('timeDelays.txt','-ascii');
nTimeDelays = length(FEDdata.timeDelays);
FEDdata.nTimeDelays = nTimeDelays;
for T = 1:nTimeDelays
    initialResiduals(T,:) = FEDdata.relativeChanges(T,:)./FEDdata.meanUncertainties;
end
FEDdata.initialResiduals = initialResiduals;
end