function [XX2,YY2] = radialDistortion(XX,YY,params)
% Input:
% XX, YY - 2d maticies of X and Y coordinates wrt distortion axis
% params - vector of distortion parameters
% aspect - ratio of distortion in X:Y
% Output:
% dXX, dYY - 2d maticies of coordinates
p = params;
R2rd = XX.^2 + YY.^2;
R = sqrt(R2rd);
XX2 = XX.*(p(1) + 0.001*p(2)*R + 1e-6*p(3)*R2rd);
YY2 = YY.*(p(1) + 0.001*p(2)*R + 1e-6*p(3)*R2rd);
end
