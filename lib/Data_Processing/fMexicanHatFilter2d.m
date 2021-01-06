function IMG2 = fMexicanHatFilter2d(IMG,sigma)

nData = size(IMG,1);
halfWindowSize = round(3*sigma);
windowSize = 2*halfWindowSize+1;
midpoint = halfWindowSize+1;
X = (1:windowSize)-midpoint;
Y = (1:windowSize)-midpoint;
[XX,YY] = meshgrid(X,Y);
RR2rd = XX.^2 + YY.^2;
sigma2rd_inv = 1/sigma^2;
ZZ = (1-0.5*RR2rd*sigma2rd_inv).*sigma2rd_inv.*exp(-RR2rd*0.5*sigma2rd_inv)/pi;
IMG2 = conv2(IMG,ZZ,'same');
end