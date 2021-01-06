function Fs = KirklandTable2f(Z,s,TableDat)
% s = 4*pi*sin(theta/2)/wavel = 2*pi*q
L0 = 3*Z - 2;
q2 = (0.5*s/pi)^2;
A = [TableDat(L0,1), TableDat(L0,3), TableDat(L0+1,1)];
B = [TableDat(L0,2), TableDat(L0,4), TableDat(L0+1,2)];
C = [TableDat(L0+1,3), TableDat(L0+2,1), TableDat(L0+2,3)];
D = [TableDat(L0+1,4), TableDat(L0+2,2), TableDat(L0+2,4)];
E = 1./(q2 + B);
term1 = A*E';
term2 = C*exp(-q2*D');
Fs = term1 + term2;
%Fs = A*E' + C*exp(-q2*D');
end