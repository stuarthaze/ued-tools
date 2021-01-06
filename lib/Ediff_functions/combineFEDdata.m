function DATA2 = combineFEDdata(DATA,indexRef)
% SAH Jan/2020
% DATAin is an array of FEDdata structures
% DATAout is a single structure of data interpolated to the temporal
% sampling in the data of indexRef

nDataSets = length(DATA);
iDat = 1:nDataSets;

% DATA2 = DATA{indexRef};
DATA2.timeDelays = DATA{indexRef}.timeDelays;
DATA2.t0 = DATA{indexRef}.t0;
DATA2.nTimeDelays = DATA{indexRef}.nTimeDelays;
xRef = DATA{indexRef}.timeDelays - DATA{indexRef}.t0;
DATA2.timeZeroCorrected = xRef;
nt = DATA2.nTimeDelays;
indx0 = 0;

% for ii = iDat(iDat~=indexRef)
for ii = iDat
    nroi = DATA{ii}.nROIs;
    xvec = DATA{ii}.timeDelays - DATA{ii}.t0;
    vroi = 1:nroi;
    
    vroi2 = vroi + indx0;
    DATA2.spotNumbers(vroi2,1) = DATA{ii}.spotNumbers;
    for roi = vroi
        % Interpolate data to new time domain
%         DATA2.I_P(:,roi+indx0) = interp1(xvec,DATA{ii}.I_P(:,roi),xRef);
%         DATA2.I_PP(:,roi+indx0) = interp1(xvec,DATA{ii}.I_PP(:,roi),xRef);
        DATA2.ratios(:,roi+indx0) = interp1(xvec,DATA{ii}.ratios(:,roi),xRef);
        DATA2.relativeChanges(:,roi+indx0) = interp1(xvec,DATA{ii}.relativeChanges(:,roi),xRef);
%         DATA2.uncertainties_P(:,roi+indx0) = interp1(xvec,DATA{ii}.uncertainties_P(:,roi),xRef);
%         DATA2.uncertainties_PP(:,roi+indx0) = interp1(xvec,DATA{ii}.uncertainties_PP(:,roi),xRef);
        DATA2.uncertainties_R(:,roi+indx0) = interp1(xvec,DATA{ii}.uncertainties_R(:,roi),xRef);
        DATA2.meanUncertainties_R(1,roi+indx0) = DATA{ii}.meanUncertainties_R(1,roi);
        
    end
    DATA2.Wts(vroi2,vroi2) = DATA{ii}.Wts;
    DATA2.HKL(vroi2,:) = DATA{ii}.HKL;
    DATA2.nROIs = indx0 + nroi;
    indx0 = DATA2.nROIs;
    % nConv - need to update later
    DATA2.relativeErrors(1,vroi2) = DATA{ii}.relativeErrors;
end
DATA2.initialResiduals = DATA2.relativeChanges./DATA2.uncertainties_R;
DATA2.rmsResid_vs_T = rms(DATA2.initialResiduals,2);
end

    