function Fit = fitTRMat_step_exp(dataMat,timeDelays)

[ny,nx] = size(dataMat);
xvec = 1:nx;
% [XX,TT] = meshgrid(1:nx,timeDelays);
Yones = ones(ny,1);
M = zeros(ny,nx);

par0(xvec) = 0;                                     % step amplitudes
par0(xvec + nx) = 0;                                % amp-exponential
par0(2*nx+1) = timeDelays(round(ny/4));             % t=0
par0(2*nx+2) = (timeDelays(end)-timeDelays(1))/5;   % sigma IRF
par0(2*nx+3) = (timeDelays(end)-timeDelays(1));     % sigma decay


fitFun = @(Pars) dataMat - params2matrix(Pars);

options = optimoptions('lsqnonlin','Display','iter','MaxFunEvals',5e4,'FunctionTolerance',1e-7);
fittedPars = lsqnonlin(fitFun,par0,[],[],options);
Fit.baseline = fittedPars(xvec);
Fit.stepAmpl = fittedPars(nx+xvec);
Fit.decayAmpl = fittedPars(2*nx+xvec);
Fit.t0 = fittedPars(3*nx+1);
Fit.sigmaIRF = fittedPars(3*nx+2);
Fit.decayConst = fittedPars(3*nx+3);
Fit.Mcalc = params2matrix(fittedPars);


    function M = params2matrix(P)
        T = timeDelays-P(2*nx+1);
        sigmaIRF = P(2*nx+2);
        sigmaExpDec = P(2*nx+3);
        unitStep = 0.5*(1+erf(T/(sigmaIRF*sqrt(2))));
        expDecay = unitStep.*(exp(-T/sigmaExpDec)-1);
        for X = 1:nx
            M(:,X) = P(X)*unitStep+P(X+nx)*expDecay;
        end
    end
end