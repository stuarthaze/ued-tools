function [bondedAtomPairs, bondedDistances] = calcBondedDistances_v3(CRYSTAL,MODEL,FitResults,cutoff)

atomNumVec = 1:size(MODEL.Fatomic,1);
atomNumbers2check = atomNumVec(MODEL.freeAtoms);
nAtFree = sum(MODEL.freeAtoms);
nDists = nAtFree*(nAtFree-1)/2;
XYZ_T = FitResults.XYZ_T;
nT = size(XYZ_T,3);

fCalc_dR = @(X) sqrt(sum(X.^2,2))';

atomPairs = zeros(nDists,2);
dists = zeros(nDists,nT);

for T = 1:nT
    XYZ = XYZ_T(:,:,T);
    L = 1;
    for indxA = 1:(nAtFree-1)
        for indxB = (indxA+1):nAtFree
            A = atomNumbers2check(indxA);
            B = atomNumbers2check(indxB);
            atomPairs(L,:) = [A, B];
            dists(L,T) = fCalc_dR(XYZ(A,:)-XYZ(B,:));
            L = L + 1;
        end
    end
end

meanDists = mean(dists,2);
bondedAtomPairsLogical = meanDists < cutoff;
bondedAtomPairs = atomPairs(bondedAtomPairsLogical,:);
bondedDistances = dists(bondedAtomPairsLogical,:);
end