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