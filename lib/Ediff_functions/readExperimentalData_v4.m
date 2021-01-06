function FED = readExperimentalData_v4(TR,INDEXING,reliabilityCutoff)

%     addpath(dirFED);
%     addpath(dirIndexing);
%     load('Wts');
%     load('HKL');
%     load('Fobs_div_Fcalc');

    indexedToUse = INDEXING.relativeErrors <= reliabilityCutoff;
    FED.spotNumbers = INDEXING.roiNumbers(indexedToUse);

    % spotIndicies_reliable = spotVec((Fobs_div_Fcalc >= reliabilityCutoffMin) & (Fobs_div_Fcalc <= reliabilityCutoffMax));
    % save('spotIndicies_reliable.mat','spotIndicies_reliable');
% 
%     relativeChanges = load('relativeChanges.txt','-ascii');
    
    FED.I_P = TR.Intensities_P(:,FED.spotNumbers);
    FED.I_PP = TR.Intensities_PP(:,FED.spotNumbers);
    FED.ratios = FED.I_PP ./FED.I_P;
    FED.relativeChanges = FED.ratios -1;
    FED.uncertainties_P = TR.Uncertainties_P(:,FED.spotNumbers);
    FED.uncertainties_PP = TR.Uncertainties_PP(:,FED.spotNumbers);
    FED.uncertainties_R = sqrt(FED.uncertainties_P.^2 + FED.uncertainties_PP.^2)./FED.I_P;
    FED.meanUncertainties_R = mean(FED.uncertainties_R,1);  
    FED.Wts = INDEXING.Wts(indexedToUse,indexedToUse); %need to re-think for num.convolutions > 1
    FED.HKL = INDEXING.HKL(indexedToUse,:);

    FED.nROIs = length(FED.spotNumbers);
    FED.nConv = 1;% INDEXING.nConv;
    FED.timeDelays = TR.timeDelays;
%     FED.Fobs_dif_Fcalc = Fobs_div_Fcalc;
%     FED.spotIndices_reliable = FED.spotNumbers;
    nTimeDelays = length(FED.timeDelays);
    FED.nTimeDelays = nTimeDelays;
    for T = 1:nTimeDelays
        initialResiduals(T,:) = FED.relativeChanges(T,:)./FED.meanUncertainties_R;
    end
    FED.initialResiduals = initialResiduals;
    FED.rmsResid_vs_T = rms(initialResiduals,2);
end