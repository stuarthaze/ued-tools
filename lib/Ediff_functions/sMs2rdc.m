function RDC = sMs2RDC(R,S,sMs,dampf)
dr = R(2)-R(1);
nR = length(R);
RDC = zeros(1,nR);
sMsDamped = sMs.*exp(-dampf*S.^2);
for r = 1:nR
    RDC(r) = sMsDamped*sin(R(r)*S');
end
end