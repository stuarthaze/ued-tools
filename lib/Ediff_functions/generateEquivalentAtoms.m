function [uvw2 atoms2] = generateEquivalentAtoms(UVW, atoms, SymStruct)
% UVW   - fractional coordinates of atoms in asymmetric unit
% atoms - atomic numbers of atoms in the asymmetric unit
    nSym = size(SymStruct.point,1);
    nAtoms = size(UVW,1);
    at_num = nAtoms;
    for sym = 1:nSym
        for atm = 1:nAtoms
            at_num = (sym-1)*nAtoms + atm;
            atoms2(at_num,1) = atoms(atm);
            uvw2(at_num,:) = UVW(atm,:).*SymStruct.point(sym,:) + SymStruct.trans(sym,:);
        end
    end
end