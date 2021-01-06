function PKS = array2peaks2013(M, minSep)

[nRows, nCols]= size(M);
PK_height_rows = zeros(nRows,nCols);
PK_height_cols = zeros(nRows,nCols);

for R = 1:nRows
    [pks, locs] = findpeaks(M(R,:),'MINPEAKDISTANCE',minSep);
    PK_height_rows(R,locs) = pks;
end

for C = 1:nCols
    [pks, locs] = findpeaks(M(:,C),'MINPEAKDISTANCE',minSep);
    PK_height_cols(locs,C) = pks;
end

PKS = PK_height_rows.*PK_height_cols;
