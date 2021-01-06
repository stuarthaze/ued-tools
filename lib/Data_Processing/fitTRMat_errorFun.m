function Fit = fitTRMat_errorFun(dataMat,timeDelays)

[ny,nx] = size(dataMat);
xvec = 1:nx;
[XX,TT] = meshgrid(1:nx,timeDelays);
Yones = ones(ny,1);

par0(1:nx) = 0;                                 % baseline
par0(nx + 1:nx) = 0;                            % step height
par0(2*nx+1) = timeDelays(round(ny/2));         % t=0
par0(2*nx+2) = timeDelays(10)-timeDelays(1);    % 

fitFun = @(Pars) dataMat - params2matrix(Pars);

fittedPars = lsqnonlin(fitFun,par0);
Fit.baseline = fittedPars(xvec);
Fit.stepHeights = fittedPars(nx+xvec);
Fit.t0 = fittedPars(2*nx+1);
Fit.sigma = fittedPars(2*nx+2);
Fit.Mcalc = params2matrix(fittedPars);


    function Mout = params2matrix(P)
        baseline = P(xvec);
        stepHeight = P(nx+xvec);
        t0 = P(2*nx+1);
        sigma = P(2*nx+2);
        baseline2d = conv2(Yones,baseline);
        stepHeights2d = conv2(Yones,stepHeight);
        unitStep = 0.5*(1+erf((TT-t0)/(sigma*sqrt(2))));
        Mout = unitStep.*stepHeights2d + baseline2d;    
    end

end