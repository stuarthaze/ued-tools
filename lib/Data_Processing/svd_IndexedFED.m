function results = svd_IndexedFED(FED)

nt = FED.nTimeDelays;
nroi = FED.nROIs;

data4svd = zeros(nroi,nt);
for roi = 1:nroi
    data4svd(roi,:) = FED.relativeChanges(:,roi)/FED.meanUncertainties_R(roi);
end

[U,S,V] = svd(data4svd);

results.U = U;
results.S = diag(S);
results.V = V;
end