function CRYST_2 = applySymToCrystal(CRYST, SymStruct)

    CRYST_2 = CRYST;
    UVW = CRYST.UVW;
    
    nSym = size(SymStruct.point,1);
    nAtoms = size(UVW,1);
    at_num = nAtoms;
    for sym = 1:nSym        
        for atm = 1:nAtoms
            at_num = (sym-1)*nAtoms + atm;
            CRYST_2.atomTypes{at_num} = CRYST.atomTypes{atm};
            CRYST_2.atomicNumbers(at_num) = CRYST.atomicNumbers(atm);
            uvw2(at_num,:) = UVW(atm,:).*SymStruct.point(sym,:) + SymStruct.trans(sym,:);
            if isfield(CRYST,'Uiso')
                CRYST_2.Uiso(at_num) = CRYST.Uiso(atm);
            end
        end
    end
    CRYST_2.UVW = uvw2;
    CRYST_2.XYZ = uvw2*CRYST.T_FracCart;
    CRYST_2.nAt = nSym*nAtoms;
end