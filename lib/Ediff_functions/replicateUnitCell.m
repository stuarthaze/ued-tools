function CRYSTAL_2 = replicateUnitCell(CRYSTAL,nRep)

CRYSTAL_2 = CRYSTAL;
UVW0 = CRYSTAL.UVW;
Uiso0 = CRYSTAL.Uiso;
nCells = nRep;
nAt = CRYSTAL.nAt;
nAt_super = nAt*nCells;
UVW = zeros(nAt_super,3);
U = zeros(nAt_super,1);
atomicNumbers = zeros(1,nAt_super);
atomNumbering = 1:nAt_super;

for cell = 1:nCells
    offset = nAt*(cell-1);
    iAtoms = (1:nAt) + offset;
    UVW(iAtoms,:) = UVW0;
    U(iAtoms) = Uiso0;
    atomicNumbers(iAtoms) = CRYSTAL.atomicNumbers;
    atomsInAsymLogical(cell,:) = (atomNumbering > offset) & (atomNumbering <= nAt*cell);

    for A=1:nAt
        atomTypes{A+offset,1} = CRYSTAL.atomTypes{A};
    end
end

CRYSTAL_2.UVW = UVW;
CRYSTAL_2.XYZ = UVW*CRYSTAL.T_FracCart;
CRYSTAL_2.Uiso = U;
CRYSTAL_2.nAt = nAt_super;
CRYSTAL_2.atomTypes = atomTypes;
CRYSTAL_2.atomicNumbers = atomicNumbers;
CRYSTAL_2.atomsInAsymLogical = atomsInAsymLogical;
CRYSTAL_2.nAt_orig = nAt;

end