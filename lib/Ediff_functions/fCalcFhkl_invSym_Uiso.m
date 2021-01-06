function Fcalc = fCalcFhkl_invSym_Uiso(ATOMS,XYZ,HKL,Fat,SYM,Uiso,modK)
% Calculate stucture factors from the asymmetric unit
% ATOMS - atomic numbers of atoms in the asymmetric unit
% XYZ   - partial atomic coordinates of atoms in the asymmetric unit
% HKL   - hkl values - (convolution number, hkl, spot number)
% Fat   - Atomic scattering factors at each contributing h,k,l, value
%         (atom number, convolution, spot)
%-------------------------------------------------------------------
[X, atomicNumbers_all, Fat_all, Uiso2] = generateEquivalentAtoms(XYZ, ATOMS, Fat, Uiso, SYM);
nConv = size(HKL,1);
nAt = size(X,1);
nSpots = size(HKL,3);
Fxyz = zeros(nConv,nAt);
XYZt = X';
Fcalc = zeros(nConv,nSpots);
k2rd = modK.^2;
Uiso2rd = Uiso.^2;
% Note: FatT(conv,atm)
% size(Fat)
% size(Fat_all)
for s = 1:nSpots
    FatT = Fat_all(:,:,s)';
    [Umesh2rd, K2rdMesh] = meshgrid(Uiso2rd,k2rd(:,s));
    T = exp(-2*pi*Umesh2rd.*K2rdMesh);   
    Fxyz(:,:) = T.*FatT.*cos(2*pi*HKL(:,:,s)*XYZt);
    Fcalc(:,s) = sum(Fxyz,2);
end
    function [xyz2, atoms2, F2, U2] = generateEquivalentAtoms(xyz, atoms, F, U, SymStruct)
        numSpots = size(F,3);
        nSym = size(SymStruct.point,1);
        nAtoms = size(xyz,1);
        xyz2 = xyz;
        at_num = 0;
        for sym = 1:nSym
            for atm = 1:nAtoms
                at_num = at_num+1;
                atoms2(at_num,1) = atoms(atm);
                xyz2(at_num,:) = xyz(atm,:).*SymStruct.point(sym,:) + SymStruct.trans(sym,:);
                for spot = 1:numSpots
                    F2(at_num,:,spot) = Fat(atm,:,spot);
                end
                U2(at_num) = U(atm);
            end
        end
    end
end
