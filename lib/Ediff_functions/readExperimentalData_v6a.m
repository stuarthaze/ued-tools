function FED = readExperimentalData_v6a(TR,INDEXING)
% Improvement on v6 - can use TR as struct or class

    FED.spotNumbers = INDEXING.roiNumbers;
    if (isa(TR,'struct') && isfield(TR,'Intensities_P')) || (~isa(TR,'struct') && ~isempty(TR.Intensities_P))
        FED.I_P = TR.Intensities_P(:,FED.spotNumbers);
    else
        FED.I_P = [];
    end
    
    if (isa(TR,'struct') && isfield(TR,'Intensities_PP')) || (~isa(TR,'struct') && ~isempty(TR.Intensities_PP))
        FED.I_PP = TR.Intensities_PP(:,FED.spotNumbers);
    else
        FED.I_PP = [];
    end
    
    if (isa(TR,'struct') && isfield(TR,'ratios')) || (~isa(TR,'struct') && ~isempty(TR.ratios))
        FED.ratios = TR.ratios(:,FED.spotNumbers);
    else
        FED.ratios = FED.I_PP ./FED.I_P;
    end
    
    FED.relativeChanges = FED.ratios -1;

    if (isa(TR,'struct') && isfield(TR,'Uncertainties_P')) || (~isa(TR,'struct') && ~isempty(TR.Uncertainties_P))
        FED.uncertainties_P = TR.Uncertainties_P(:,FED.spotNumbers);
    else
        FED.uncertainties_P = [];
    end
    
    if (isa(TR,'struct') && isfield(TR,'Uncertainties_PP')) || (~isa(TR,'struct') && ~isempty(TR.Uncertainties_PP))    
        FED.uncertainties_PP = TR.Uncertainties_PP(:,FED.spotNumbers);
    else
        FED.uncertainties_PP = [];
    end
    
    if (isa(TR,'struct') && isfield(TR,'Uncertainties_R')) || (~isa(TR,'struct') && ~isempty(TR.Uncertainties_R))    
        FED.uncertainties_R = TR.Uncertainties_R(:,FED.spotNumbers);
    else
        FED.uncertainties_R = sqrt(FED.uncertainties_P.^2 + FED.uncertainties_PP.^2)./FED.I_P;
    end
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
    FED.relativeErrors = INDEXING.relativeErrors';

end