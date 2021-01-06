function CRYSTAL_2 = createSuperCell(CRYSTAL,supercellDimensions)

CRYSTAL_2 = CRYSTAL;
XYZ0 = CRYSTAL.XYZ;
Uiso0 = CRYSTAL.Uiso;
nx = supercellDimensions(1);
ny = supercellDimensions(2);
nz = supercellDimensions(3);

nCells = nx*ny*nz;
nAt = CRYSTAL.nAt;
nAt_super = nAt*nCells;
XYZ = zeros(nAt_super,3);
U = zeros(nAt_super,1);
atomicNumbers = zeros(1,nAt_super);

iter = 0;
for x = 1:nx
    for y = 1:ny
        for z = 1:nz
            iter = iter + 1;
            offset = nAt*(iter-1);
            iAtoms = (1:nAt) + offset;
            Tvec = [x-1,y-1,z-1];
            [T2d,Indx] = meshgrid(Tvec,1:nAt);
            XYZ(iAtoms,:) = XYZ0 + T2d;
            U(iAtoms) = Uiso0;
            atomicNumbers(iAtoms) = CRYSTAL.atomicNumbers;
            for A=1:nAt
                atomTypes{A+offset,1} = CRYSTAL.atomTypes{A};
            end
        end
    end
end


CRYSTAL_2.XYZ = XYZ;
CRYSTAL_2.Uiso = U;
CRYSTAL_2.nAt = nAt_super;
CRYSTAL_2.atomTypes = atomTypes;
CRYSTAL_2.atomicNumbers = atomicNumbers;

end