function F = calculateFrequenciesFFT(dt,n)

fs = 1/dt;
df = fs/n;
F0 = ((0:n-1)-floor(n/2))*df;
F = ifftshift(F0);
end