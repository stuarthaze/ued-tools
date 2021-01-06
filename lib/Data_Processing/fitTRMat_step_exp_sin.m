function Fit = fitTRMat_step_exp_sin(dataMat,timeDelays)

[ny,nx] = size(dataMat);
xvec = 1:nx;
% [XX,TT] = meshgrid(1:nx,timeDelays);
Yones = ones(ny,1);
dt = timeDelays(2)-timeDelays(1);

par0(xvec) = 0;                                     % step amplitudes
par0(xvec + nx) = 0;                                % amp-exponential
par0(xvec + 2*nx) = 0;                              % amp sin
par0(3*nx+1) = timeDelays(round(ny/3));             % t=0
par0(3*nx+2) = (timeDelays(end)-timeDelays(1))/10;  % sigma IRF
par0(3*nx+3) = (timeDelays(end)-timeDelays(1));     % sigma decay
par0(3*nx+4) = 0.6e-3;                              % oscillation freq in fs-1
par0(3*nx+5) = (timeDelays(end)-timeDelays(1));     % oscillation decay


fitFun1 = @(Pars) dataMat - params2matrix1(Pars);
fitFun2 = @(Pars) dataMat - params2matrix2(Pars);

options = optimoptions('lsqnonlin','Display','iter','MaxFunEvals',5e4,'FunctionTolerance',1e-7);
fittedPars1 = lsqnonlin(fitFun1,par0,[],[],options);
disp('fitted step, now add additional functions');
fittedPars2 = lsqnonlin(fitFun2,fittedPars1,[],[],options);
Fit.baseline = fittedPars2(xvec);
Fit.stepAmpl = fittedPars2(nx+xvec);
Fit.decayAmpl = fittedPars2(2*nx+xvec);
Fit.t0 = fittedPars2(3*nx+1);
Fit.sigmaIRF = fittedPars2(3*nx+2);
Fit.expDecayConst = fittedPars2(3*nx+3);
Fit.sinFreq = fittedPars2(3*nx+4);
Fit.sinDecayConst = fittedPars2(3*nx+5);
Fit.Mcalc = params2matrix2(fittedPars2);


    function M = params2matrix1(P)
        M = zeros(ny,nx);
        T = timeDelays-P(3*nx+1);
        sigmaIRF = P(3*nx+2);
        sigmaExpDec = P(3*nx+3);
        f = P(3*nx+4);
        sinDecayConst = P(3*nx+5);
        unitStep = 0.5*(1+erf(T/(sigmaIRF*sqrt(2))));
        expDecay = unitStep.*(exp(-T/sigmaExpDec)-1);
        unitSin = unitStep.*exp(-T/sinDecayConst).*sin(2*pi*f*T);
        for X = 1:nx
            M(:,X) = P(X)*unitStep;
        end
    end

    function M = params2matrix2(P)
        M = zeros(ny,nx);
        T = timeDelays-P(3*nx+1);
        Tpos = T >= 0;
        sigmaIRF = P(3*nx+2);
        expDecayConst = P(3*nx+3);
        f = P(3*nx+4);
        sinDecayConst = P(3*nx+5);
        unitStep = 0.5*(1+erf(T/(sigmaIRF*sqrt(2))));
        expDecay = Tpos.*(exp(-T/expDecayConst)-1);
        unitSin = Tpos.*exp(-T/sinDecayConst).*sin(2*pi*f*T);
        for X = 1:nx
            M(:,X) = P(X)*unitStep + fGaussianFilter(P(X+2*nx)*unitSin + P(X+nx)*expDecay, sigmaIRF/dt);
        end
    end

    function F = fGaussianFilter(DATA,sigma)
        nDat = length(DATA);
        halfWindowSize = round(3*sigma);
        windowSize = 2*halfWindowSize+1;
        midpoint = halfWindowSize+1;
        X = (1:windowSize)-midpoint;
        Y = exp(-0.5*X.^2./sigma^2)/(sqrt(2*pi)*sigma);
        DATApad(1:halfWindowSize) = DATA(end);
        DATA2 = [DATA,DATApad];
        Fpad = conv(DATA2,Y,'same');
        F = Fpad(1:nDat);
    end

end