function Iat = calcIat(Fs)
nat = length(Fs(:,1));
for a = 1:nat
    IatA(a,:) = Fs(a,:).^2;
end
Iat = sum(IatA,1);
end
