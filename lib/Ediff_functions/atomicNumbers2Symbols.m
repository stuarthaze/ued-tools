function atomSymbols = atomicNumbers2Symbols(atomicNumbers)
nAtoms = length(atomicNumbers);
element{1} = 'H';
element{2} = 'He';
element{3} = 'Li';
element{4} = 'Be';
element{6} = 'C';
element{15} = 'P';
element{16} = 'S';
element{34} = 'Se';
element{50} = 'Sn';
element{78} = 'Pt';
for atm = 1:nAtoms
    atomSymbols{atm} = element(atomicNumbers(atm));
end
end