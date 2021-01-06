function Fs = KirklandTable2f_v2(Z,s,table)
% v2 - vectorized
L0 = 3*Z - 2;
q2 = (0.5*s/pi).^2;
A = [table(L0,1) table(L0,3) table(L0+1,1)];
B = [table(L0,2) table(L0,4) table(L0+1,2)];
C = [table(L0+1,3) table(L0+2,1) table(L0+2,3)];
D = [table(L0+1,4) table(L0+2,2) table(L0+2,4)];
E = 1./(q2 + B);
term1 = A*E';
term2 = C*exp(-q2*D');
Fs = term1 + term2;
%Fs = A*E' + C*exp(-q2*D');
end