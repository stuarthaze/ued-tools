function Fcalc = fXYZ2Fspots(XYZ,HKL,Fat)
% Calculates F_at for a set of atoms for given HKL values
% Fat - atomic structure factor - dimensions [nAt, nConv, nSpots]
nConv = size(HKL,1);
nAt = size(XYZ,1);
nSpots = size(HKL,3);
Fxyz = zeros(nConv,nAt);
XYZt = XYZ';
Fcalc = zeros(nConv,nSpots);
for s = 1:nSpots
    FatT = Fat(:,:,s)';
    Fxyz(:,:) = FatT.*exp(2*pi*1i*HKL(:,:,s)*XYZt);
    Fcalc(:,s) = sum(Fxyz,2);
end
end