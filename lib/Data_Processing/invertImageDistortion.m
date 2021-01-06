function [Xout, Yout] = invertImageDistortion(Xdist,Ydist)
% SHayes 2020
% Input: Xdist,Ydist - 2d matricies for distortion
% Output: Xout, Yout - 2d matricies for reversing distortion

[ny,nx] = size(Xdist);
[Xmesh,Ymesh] = meshgrid(1:nx,1:ny);
xstart = [Xmesh,Ymesh];
X1 = interp2(Xmesh,Xdist,Ydist);
Y1 = interp2(Ymesh,Xdist,Ydist);
options = optimset('Display','iter','MaxFunEvals',1e6);
% fun2fit = @(P) sum(sum(calcErrors(P).^2));
% xfit = fminunc(fun2fit,xstart,options);
xfit = lsqnonlin(@calcErrors,xstart,[],[],options);
Xout = xfit(:,1:nx);
Yout = xfit(:,nx+(1:nx));

    function errorMat = calcErrors(M)
        Xerr = interp2(X1,M(:,1:nx),M(:,(1:nx)+nx))-Xmesh;
        Yerr = interp2(Y1,M(:,1:nx),M(:,(1:nx)+nx))-Ymesh;
        
        errorMat = [Xerr,Yerr];
        errorMat(isnan(errorMat)) = 0;
    end
end

