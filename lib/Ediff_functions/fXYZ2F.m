function Fcalc = fXYZ2F(XYZ,HKL,Fat)
% Calculates F at for a set of atoms for given HKL values
% Fat - atomic structure factor - dimensions [nAt x nHKL]
nHKL = size(HKL,1);
FatT = Fat';
nAt = size(XYZ,1);
Fxyz = zeros(nHKL,nAt);
HKLt = HKL';
XYZt = XYZ';
Fxyz = FatT.*exp(2*pi*i*HKL*XYZt);
Fcalc = sum(Fxyz,2);