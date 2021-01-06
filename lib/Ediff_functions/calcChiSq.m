function chiSq = calcChiSq(sMsExp,sMsTheor,sMsExpSigmas)
   ns = length(sMsExp);
   vecChiSq = (sMsExp-sMsTheor).^2./(sMsExpSigmas.^2);
   chiSq = sum(vecChiSq)/ns;
end