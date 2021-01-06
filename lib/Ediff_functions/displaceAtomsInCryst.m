function CRYST2 = displaceAtomsInCryst(CRYST,dispUVW,atomNumbers)

UVW2 = CRYST.UVW;
for ii = 1:length(atomNumbers)
    atm = atomNumbers(ii);
    UVW2(atm,:) = CRYST.UVW(atm,:) + dispUVW;
end
CRYST2 = CRYST;
CRYST2.UVW = UVW2;
CRYST2.XYZ = UVW2*CRYST.T_FracCart;
end
