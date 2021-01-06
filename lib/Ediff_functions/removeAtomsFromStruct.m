function cryst2 = removeAtomsFromStruct(cryst,atoms2remove)

natoms = cryst.nAt;

if islogical(atoms2remove)
    atoms2keep = ~atoms2remove;
else
    nDup = length(atoms2remove);
    atoms2keep = logical(1:natoms);
    for ii = 1:nDup
        atoms2keep(atoms2remove(ii)) = 0;
    end
end

cryst2 = cryst;
cryst2.UVW = cryst.UVW(atoms2keep,:);
cryst2.atomTypes = cryst.atomTypes(atoms2keep);
cryst2.atomicNumbers = cryst.atomicNumbers(atoms2keep);
cryst2.XYZ = cryst.XYZ(atoms2keep,:);
cryst2.nAt = sum(atoms2keep);
cryst2.Uiso = cryst.Uiso(atoms2keep);

end
