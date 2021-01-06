function IMG2 = fGaussianFilter2d(IMG,sigma)

halfWindowSize = round(3*sigma);
windowSize = 2*halfWindowSize+1;
midpoint = halfWindowSize+1;
X = (1:windowSize)-midpoint;
Y = (1:windowSize)-midpoint;
[XX,YY] = meshgrid(X,Y);
RR2rd = XX.^2 + YY.^2;
ZZ = exp(-RR2rd*0.5/sigma^2);
A = sum(sum(ZZ));
ZZscaled = ZZ/A;

IMG2 = conv2(IMG,ZZscaled,'same');
end