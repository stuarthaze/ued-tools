function Radii = readCovalentRadiiFromFile(atomicNumbers)

    data = xlsread('CovalentRadii_Codero2008DaltonTrans.xlsx');

    columnRadii = 3;
    columnAtomicNumber = 1;
    Radii = zeros(size(atomicNumbers));
    nAt = length(atomicNumbers);
    for atm = 1:nAt
        Radii(atm) = data(atomicNumbers(atm)==data(:,columnAtomicNumber),columnRadii);
    end
end