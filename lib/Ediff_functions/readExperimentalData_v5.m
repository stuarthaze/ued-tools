function FED = readExperimentalData_v5(TR,INDEXING)

    FED.spotNumbers = INDEXING.roiNumbers;
    
    FED.I_P = TR.Intensities_P(:,FED.spotNumbers);
    FED.I_PP = TR.Intensities_PP(:,FED.spotNumbers);
    FED.ratios = FED.I_PP ./FED.I_P;
    FED.relativeChanges = FED.ratios -1;
    FED.uncertainties_P = TR.Uncertainties_P(:,FED.spotNumbers);
    FED.uncertainties_PP = TR.Uncertainties_PP(:,FED.spotNumbers);
    FED.uncertainties_R = sqrt(FED.uncertainties_P.^2 + FED.uncertainties_PP.^2)./FED.I_P;
    FED.meanUncertainties_R = mean(FED.uncertainties_R,1);  
    FED.Wts = INDEXING.Wts;
    FED.HKL = INDEXING.HKL;

    FED.nROIs = length(FED.spotNumbers);
    FED.nConv = 1;% INDEXING.nConv;
    FED.timeDelays = TR.timeDelays;
    nTimeDelays = length(FED.timeDelays);
    FED.nTimeDelays = nTimeDelays;
    for T = 1:nTimeDelays
        initialResiduals(T,:) = FED.relativeChanges(T,:)./FED.meanUncertainties_R;
    end
    FED.initialResiduals = initialResiduals;
    FED.rmsResid_vs_T = rms(initialResiduals,2);

end