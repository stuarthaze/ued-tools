function chiSq = calcChiSq(sMsExp,sMsTheor)
   ns = length(sMsExp);
   vecChiSq = (sMsExp-sMsTheor).^2;
   chiSq = sum(vecChiSq)/ns;
end