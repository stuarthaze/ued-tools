function [PKS_h, PKS_w, PKS_p] = array2peaks2017(M)

[nRows, nCols]= size(M);
PK_height_rows = zeros(nRows,nCols);
PK_width_rows = zeros(nRows,nCols);
PK_prominance_rows = zeros(nRows,nCols);

PKS_height_cols = zeros(nRows,nCols);
PKS_width_cols = zeros(nRows,nCols);
PKS_prominance_cols = zeros(nRows,nCols);

for R = 1:nRows
    [pks, locs, wd, pr] = findpeaks(M(R,:));
    PK_height_rows(R,locs) = pks;
    PK_width_rows(R,locs) = wd;
    PK_prominance_rows(R,locs) = pr;
end

for C = 1:nCols
    [pks, locs, wd, pr] = findpeaks(M(:,C));
    PK_height_cols(locs,C) = pks;
    PK_width_cols(locs,C) = wd;
    PK_prominance_cols(locs,C) = pr;
end

PKS_h = PK_height_rows + PK_height_cols;
PKS_w = PK_width_rows + PK_width_cols;
PKS_p = PK_prominance_rows + PK_prominance_cols;