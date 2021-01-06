function amplitudes = generateApproxAmplitudes(distList)
% Generate approximate vibrational amplitudes according to distance
% Amp = 0.05 + 0.025*r
amplitudes = 0.05 + 0.025*distList(:,3);
end