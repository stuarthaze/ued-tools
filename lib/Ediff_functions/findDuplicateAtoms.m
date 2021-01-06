function duplicates = findDuplicateAtoms(CRYST,rTolerance)

cutoff = rTolerance^2;
nAt = size(CRYST.UVW,1);
xyz = CRYST.UVW*CRYST.T_FracCart;
Duplicates2d = zeros(nAt);
for A = 1:(nAt-1)
    for B = (A+1):nAt
        Rab_2rd = sum((xyz(A,:)-xyz(B,:)).^2);
        Duplicates2d(A,B) = Rab_2rd <= cutoff;
    end
end
duplicates = logical(sum(Duplicates2d,1));

end