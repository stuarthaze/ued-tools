function [xyz2 atoms2 occupancies2] = generateEquivalentAtomsOccupancies(xyz, atoms, SymStruct,occupancies)
% xyz   - fractional coordinates of atoms in asymmetric unit
% atoms - atomic numbers of atoms in the asymmetric unit
    nSym = size(SymStruct.point,1);
    nAtoms = size(xyz,1);
    at_num = nAtoms;
    for sym = 1:nSym
        for atm = 1:nAtoms
            at_num = (sym-1)*nAtoms + atm;
            atoms2(at_num,1) = atoms(atm);
            xyz2(at_num,:) = xyz(atm,:).*SymStruct.point(sym,:) + SymStruct.trans(sym,:);
            occupancies2(at_num) = occupancies(atm);
        end
    end
end