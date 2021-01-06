function DATA2 = replaceNaNsWithNearestValue(DATA)
nData = length(DATA);
X = 1:nData;
DATA2 = DATA;

NANs = isnan(DATA);
NANindicies = X(NANs);
notNaNs = ~NANs;
notNaNindicies = X(notNaNs);
numNaNs = sum(NANs);

if numNaNs
    for y = 1:numNaNs
        IndxDiff = abs(NANindicies(y)-notNaNindicies);
        [value, indx] = min(IndxDiff);
        DATA2(NANindicies(y)) = DATA(notNaNindicies(indx));
    end
end
end