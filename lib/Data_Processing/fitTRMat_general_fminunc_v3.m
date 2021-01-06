function FitOut = fitTRMat_general_fminunc_v3(dataMat,timeDelays,fitoptions)
% fitoptions includes:
%   numExp
%   numSin - if not==0 also fitoptions.freqStart


numExp = fitoptions.numExp;
numSin = fitoptions.numSin;
if numSin
    freqStart = fitoptions.freqStart;
    if numSin ~= length(freqStart)
        disp('Number of sin components is different from the number of input frequencies');
    end
end

[ny,nx] = size(dataMat);
xvec = 1:nx;
dt = timeDelays(2)-timeDelays(1);

% ---- Calculate Parameter Indices for each term -----------
ParIndx_Stp = xvec;
for indxExp = 1:numExp
    ParIndx_Exp(indxExp,xvec) = indxExp*nx + xvec;
end
for indxSin = 1:numSin
    ParIndx_Sin(indxSin,xvec) = numExp*nx + indxSin*nx + xvec;
end
ParIndx_T0 = nx + numExp*nx + numSin*nx + 1;
ParIndx_IRF = ParIndx_T0 + 1;
ParIndx_ExpConst = ParIndx_IRF + (1:numExp);
ParIndx_Freq = ParIndx_IRF + numExp + (1:numSin);
ParIndx_SinDecays = ParIndx_IRF + numExp + numSin + (1:numSin);
%--------------------------------------------------------------
% Set starting values
Val0(ParIndx_Stp) = 0;                                          % step amplitudes
if numExp
    Val0(ParIndx_Exp) = 0;                                      % amp exponential
    Val0(ParIndx_ExpConst) = (timeDelays(end)-timeDelays(1));   % exponential decay
end
if numSin
    Val0(ParIndx_Sin) = 0;                                      % amp sin
    Val0(ParIndx_Freq) = freqStart;
    Val0(ParIndx_SinDecays) = (timeDelays(end)-timeDelays(1))/2;
end    
Val0(ParIndx_T0) = timeDelays(round(ny/4));             
Val0(ParIndx_IRF) = (timeDelays(end)-timeDelays(1))/10; 

% Use to display starting values:
StartVals = params2struct(Val0)
%-------------------------------------------------
% Fitting
fitFun = @(Pars) sum(rms(dataMat - params2matrix_gen(Pars)));

options = optimoptions('fminunc','Algorithm','quasi-newton','Display','iter','MaxFunEvals',1e5,'FunctionTolerance',1e-7);
% First fit only step, t0 and IRF
numExp_in = numExp;
numSin_in = numSin;
numExp = 0;
numSin = 0;
fittedPars1 = fminunc(fitFun,Val0,options);
% Display results of initial fit
Fit1 = params2struct(fittedPars1)
% Add sin and exp decay components
numExp = numExp_in;
numSin = numSin_in;
fittedPars2 = fminunc(fitFun,fittedPars1,options);
Fit = params2struct(fittedPars2);
%----------------------------------
% OUTPUT
FitOut = Fit;
FitOut.numExp = fitoptions.numExp;
FitOut.numSin = fitoptions.numSin;
FitOut.Mdata = dataMat;
FitOut.Mcalc = params2matrix_gen(fittedPars2);
FitOut.Mdiff = dataMat - FitOut.Mcalc;
FitOut.timeDelays = timeDelays;
%------------------------------------
    function M = params2matrix_gen(P)
        M = zeros(ny,nx);
        decayMat = zeros(ny,nx);
        sinMat = zeros(ny,nx);
        T = timeDelays-P(ParIndx_T0);
        Tpos = T >= 0;
        for ii = 1:numExp
            P(ParIndx_ExpConst(ii));
%             unitDecay =Tpos.*(exp(-T/P(ParIndx_ExpConst(ii)))-1);
            unitDecay = imgaussfilt(Tpos.*(exp(-T/P(ParIndx_ExpConst(ii)))-1),P(ParIndx_IRF)/dt);
            [decayAmps2d, unitDecay2d] = meshgrid(P(ParIndx_Exp(ii,:)),unitDecay);
            decayMat = decayMat + decayAmps2d.*unitDecay2d;
        end
        for ii = 1:numSin
            unitSin = imgaussfilt(Tpos.*(exp(-T/P(ParIndx_SinDecays(ii))).*sin(2*pi*P(ParIndx_Freq(ii))*T)),P(ParIndx_IRF)/dt);
            [sinAmps2d, unitSin2d] = meshgrid(P(ParIndx_Sin(ii,:)),unitSin);
            sinMat = sinMat + sinAmps2d.*unitSin2d;
        end
        unitStep = 0.5*(1+erf(T/(P(ParIndx_IRF)*sqrt(2))));
        [stepAmps2d,unitStep2d] = meshgrid(P(ParIndx_Stp),unitStep);
        if numExp && numSin
            M = stepAmps2d.*unitStep2d + decayMat + sinMat;
        elseif numExp
            M = stepAmps2d.*unitStep2d + decayMat;
        elseif numSin
            M = stepAmps2d.*unitStep2d + sinMat;
        else
            M = stepAmps2d.*unitStep2d;
        end
    end

    function STRC = params2struct(P)
        STRC.stepAmpl = P(ParIndx_Stp);
        if numExp
            STRC.decayAmpl = P(ParIndx_Exp);
            STRC.expDecayConst = P(ParIndx_ExpConst);
        end
        if numSin
            STRC.sinAmpl = P(ParIndx_Sin);
            STRC.sinFreq = P(ParIndx_Freq);
            STRC.sinDecayConst = P(ParIndx_SinDecays);
        end
        STRC.t0 = P(ParIndx_T0);
        STRC.sigmaIRF = P(ParIndx_IRF);
        STRC.Mcalc = params2matrix_gen(P);
    end

end