function CRYSTAL_2 = createSuperCell_v3(CRYSTAL,supercellDimensions)

CRYSTAL_2 = CRYSTAL;
UVW0 = CRYSTAL.UVW;
Uiso0 = CRYSTAL.Uiso;
na = supercellDimensions(1);
nb = supercellDimensions(2);
nc = supercellDimensions(3);

nCells = na*nb*nc;
nAt = CRYSTAL.nAt;
nAt_super = nAt*nCells;
UVW = zeros(nAt_super,3);
U = zeros(nAt_super,1);
atomicNumbers = zeros(1,nAt_super);

iter = 0;
for x = 1:na
    for y = 1:nb
        for z = 1:nc
            iter = iter + 1;
            offset = nAt*(iter-1);
            iAtoms = (1:nAt) + offset;
            Tvec = [x-1,y-1,z-1];
            [T2d,Indx] = meshgrid(Tvec,1:nAt);
            UVW(iAtoms,:) = UVW0 + T2d;
            U(iAtoms) = Uiso0;
            atomicNumbers(iAtoms) = CRYSTAL.atomicNumbers;
            for A=1:nAt
                atomTypes{A+offset,1} = CRYSTAL.atomTypes{A};
            end
        end
    end
end


CRYSTAL_2.UVW = UVW;
CRYSTAL_2.XYZ = UVW*CRYSTAL.T_FracCart;
CRYSTAL_2.Uiso = U;
CRYSTAL_2.nAt = nAt_super;
CRYSTAL_2.atomTypes = atomTypes;
CRYSTAL_2.atomicNumbers = atomicNumbers;

end