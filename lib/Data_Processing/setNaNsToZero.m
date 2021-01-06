function M2 = setNaNsToZero(M)

M2 = M;
M2(isnan(M)) = 0;
end