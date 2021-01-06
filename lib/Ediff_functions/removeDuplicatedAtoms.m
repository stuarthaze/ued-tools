function cryst2 = removeDuplicatedAtoms(cryst,duplicates)

cryst2 = cryst;
indx2keep = ~duplicates;
cryst2.UVW = cryst.UVW(indx2keep,:);
cryst2.XYZ = cryst.XYZ(indx2keep,:);
cryst2.nAt = sum(indx2keep);
cryst2.atomicNumbers = cryst.atomicNumbers(indx2keep);
cryst2.atomTypes = [];
if length(cryst.Uiso) == cryst.nAt
    cryst2.Uiso = cryst.Uiso(indx2keep);
end
end