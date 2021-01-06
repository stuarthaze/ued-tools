function F_at = calculateAtomicF_Kirkland_v2(CrystalStructure,HKLs)

ScatFactTable = load('KirklandTable','-ascii');
nHKL = size(HKLs,1);
nAt = CrystalStructure.nAt;
% Atomic Scattering Factors
modS = 2*pi*fHKL2modK(HKLs,CrystalStructure.axes);
F_at = zeros([nAt,nHKL]);
for s = 1:nHKL    
    for atm = 1:nAt
        F_at(atm,s) = KirklandTable2f(CrystalStructure.atomicNumbers(atm),modS(s),ScatFactTable);
    end
end
end