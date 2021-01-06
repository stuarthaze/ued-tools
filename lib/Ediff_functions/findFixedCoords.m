function fixedCoords = findFixedCoords(XYZ,SYMstruct)
% Not working properly - need to use symmetry operations for special
% rather than general coordinates.
nSym = size(SYMstruct.point,1);
nAt = size(XYZ,1);
symPoints = SYMstruct.trans./(1-SYMstruct.point);
for N = 1:nAt
    for sym = 1:nSym
        atomFixed(sym,:) = XYZ(N,:) == symPoints(sym,:);
    end
    fixedCoords(N,:) = max(atomFixed,[],1);
end
end