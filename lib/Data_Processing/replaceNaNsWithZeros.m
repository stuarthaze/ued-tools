function DATA2 = replaceNaNsWithZeros(DATA)
    DATA2 = DATA;
    NANs = isnan(DATA);
    DATA2(NANs) = 0;
end