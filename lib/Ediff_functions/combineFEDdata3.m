function DATA2 = combineFEDdata3(DATA,Tnew,data2use)
% SAH Jan/2020
% DATAin is an array of FEDdata structures
% DATAout is a single structure of data interpolated to the temporal
% tVec is the new sampling for interpolation
% [data2use] optional - define indices of data to combine
% v3 - Weighted uncertainties according to time step (dt)

nDataSets = length(DATA);
if nargin == 3
    iDat = data2use;
else
    iDat = 1:nDataSets;
end

% DATA2 = DATA{indeTnew};
DATA2.timeDelays = Tnew;
DATA2.t0 = 0;
DATA2.nTimeDelays = length(Tnew);
DATA2.timeZeroCorrected = Tnew;
indx0 = 0;

% Calculate time steps
dt_new = (Tnew(end)-Tnew(1))/(length(Tnew)-1);

% 
% for ii = iDat(iDat~=indeTnew)
for ii = iDat
    nroi = DATA{ii}.nROIs;
    xvec = DATA{ii}.timeDelays - DATA{ii}.t0;
    dt = (xvec(end)-xvec(1))/(length(xvec)-1);
    
    errorScaling = sqrt(dt/dt_new);
    if errorScaling < 1
        errorScaling = 1;
    end
    
    vroi = 1:nroi;
    vroi2 = vroi + indx0;
    DATA2.spotNumbers(vroi2,1) = DATA{ii}.spotNumbers;
    
    for roi = vroi
        % Interpolate data to new time domain
%         DATA2.I_P(:,roi+indx0) = interp1(xvec,DATA{ii}.I_P(:,roi),Tnew);
%         DATA2.I_PP(:,roi+indx0) = interp1(xvec,DATA{ii}.I_PP(:,roi),Tnew);
        DATA2.ratios(:,roi+indx0) = interp1(xvec,DATA{ii}.ratios(:,roi),Tnew);
        DATA2.relativeChanges(:,roi+indx0) = interp1(xvec,DATA{ii}.relativeChanges(:,roi),Tnew);
%         DATA2.uncertainties_P(:,roi+indx0) = interp1(xvec,DATA{ii}.uncertainties_P(:,roi),Tnew);
%         DATA2.uncertainties_PP(:,roi+indx0) = interp1(xvec,DATA{ii}.uncertainties_PP(:,roi),Tnew);
        DATA2.uncertainties_R(:,roi+indx0) = interp1(xvec,DATA{ii}.uncertainties_R(:,roi),Tnew)*errorScaling;
        DATA2.meanUncertainties_R(1,roi+indx0) = DATA{ii}.meanUncertainties_R(1,roi)*errorScaling;
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

    