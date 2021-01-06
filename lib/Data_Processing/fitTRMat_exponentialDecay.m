function Fit = fitTRMat_exponentialDecay(dataMat,timeDelays)

[ny,nx] = size(dataMat);
xvec = 1:nx;
[XX,TT] = meshgrid(1:nx,timeDelays);
Yones = ones(ny,1);
M = zeros(ny,nx);

par0(1:nx) = 0;                                     % baseline
par0(nx + 1:nx) = 0;                                % amplitudes
par0(2*nx+1) = timeDelays(round(ny/2));             % t=0
par0(2*nx+2) = (timeDelays(end)-timeDelays(1))/5;   % sigma 

fitFun = @(Pars) dataMat - params2matrix2(Pars);

options = optimoptions('lsqnonlin','Display','iter','MaxFunEvals',5e4,'FunctionTolerance',1e-7);
fittedPars = lsqnonlin(fitFun,par0,[],[],options);
Fit.baseline = fittedPars(xvec);
Fit.amplitude = fittedPars(nx+xvec);
Fit.t0 = fittedPars(2*nx+1);
Fit.sigma = fittedPars(2*nx+2);
Fit.Mcalc = params2matrix2(fittedPars);

% 
%     function Mout = params2matrix(P)
%         baseline = P(xvec);
%         amplitudes = P(nx+xvec);
%         t0 = P(2*nx+1);
%         sigma = P(2*nx+2);
%         baseline2d = conv2(Yones,baseline);
%         amplitudes2d = conv2(Yones,amplitudes);
%         TTpos = TT >= t0;
%         unitExponential = 1-exp(-(TT-t0)/sigma);
%         Mout = TTpos.*unitExponential.*amplitudes2d + baseline2d;    
%     end

    function M = params2matrix2(P)
        T = timeDelays-P(2*nx+1);
        Tpos = T > 0;
        expDecay = Tpos.*(1-exp(-T/P(2*nx+2)));
        for X = 1:nx
            M(:,X) = P(X)+P(X+nx)*expDecay;
        end
    end
end