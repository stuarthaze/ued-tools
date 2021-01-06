function [vecMajor, vecMinor] = calcMajorMinorAxes(vecNorm,vecApproxMaj)
dirMinor = cross(vecNorm,vecApproxMaj);
vecMinor = dirMinor/sqrt(sum(dirMinor.^2));
dirMajor = cross(vecMinor,vecNorm);
vecMajor = dirMajor/sqrt(sum(dirMajor.^2));
end
