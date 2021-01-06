function Iout = calcPowderDiffraction1d(qvec,qrings,intensities,width)

    nRings = length(qrings);
    nq = length(qvec);
    Imat = zeros([nRings,nq]);

    for ii = 1:nRings
        ds = qvec - qrings(ii);
        Imat(ii,:) = intensities(ii)*exp(-0.5*ds.^2/(width^2));
    end
    Iout = sum(Imat,1);
end