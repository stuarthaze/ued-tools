function F_at = calculateAtomicF_Kirkland(CrystalStructure,HKLs)

    ScatFactTable = load('KirklandTable','-ascii');
    nConv = size(HKLs,1);
    nSpots = size(HKLs,3);
    nAt = CrystalStructure.nAt;
    % Atomic Scattering Factors
    for s = 1:nSpots
        modS(:,s) = 2*pi*fHKL2modK(HKLs(:,:,s),CrystalStructure.axes);
        for atm = 1:nAt
            for conv = 1:nConv
                F_at(atm,conv,s) = KirklandTable2f(CrystalStructure.atomicNumbers(atm),modS(conv,s),ScatFactTable);
            end
        end
    end
    function mod_k = fHKL2modK(HKL,axes)
        nSpt = size(HKL,3);
        nConv = size(HKL,1);
        mod_k = zeros(nConv,nSpt);
        abc_recip = zeros(3);
        abc_recip(1,:) = axes.a_star;
        abc_recip(2,:) = axes.b_star;
        abc_recip(3,:) = axes.c_star;
        for S = 1:nSpt
            R_star = HKL(:,:,S)*abc_recip;
            mod_k(:,S) = sqrt(sum(R_star.^2,2));
        end
    end
end