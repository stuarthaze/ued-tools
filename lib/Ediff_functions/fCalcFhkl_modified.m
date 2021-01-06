function Fcalc = fCalcFhkl(ATOMS,XYZ,HKL,Fat,SYM)
% Calculate stucture factors from the asymmetric unit
% ATOMS - atomic numbers of atoms in the asymmetric unit
% XYZ   - partial atomic coordinates of atoms in the asymmetric unit
% HKL   - hkl values - (convolution number, hkl, spot number)
% Fat   - Atomic scattering factors at each contributing h,k,l, value
%         (atom number, convolution, spot)
%-------------------------------------------------------------------
[X atoms Fat_all] = generateEquivalentAtoms(XYZ, ATOMS, SYM);
nConv = size(HKL,1);
nAt = size(X,1);
nSpots = size(HKL,3);
Fxyz = zeros(nConv,nAt);
XYZt = X';
Fcalc = zeros(nConv,nSpots);
for s = 1:nSpots
    FatT = Fat_all(:,:,s)';
    Fxyz(:,:) = FatT.*exp(2*pi*1i*HKL(:,:,s)*XYZt);
    Fcalc(:,s) = sum(Fxyz,2);
end
    function [xyz2 atoms2 F] = generateEquivalentAtoms(xyz, atoms, SymStruct)
        nSym = size(SymStruct.point,1);
        symElements = SymStruct.trans./(1-SymStruct.point);
        indices2correct = find((SymStruct.trans==0)&(SymStruct.point==1));
        symElements(indices2correct) = 1;
        nAtoms = size(xyz,1);
        xyz2 = xyz;
        at_num = nAtoms;
        for sym = 2:nSym
            for atm = 1:nAtoms
                if (xyz(atm,:) ~= symElements(sym,:))|(symElements(sym,:)==1)
                    at_num = at_num+1;
                    atoms2(at_num,1) = atoms(atm);
                    xyz2(at_num,:) = xyz(atm,:).*SymStruct.point(sym,:) + SymStruct.trans(sym,:);
                    F(at_num,:,:) = Fat(atm,:,:);
                end
            end
        end
    end
end
