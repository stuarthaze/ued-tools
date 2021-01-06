function Ispot = fCalcSpotIntensities(FcalcSpots,Wts)
    Ispot = sum(FcalcSpots.^2.*Wts, 1);
end