function F = fGaussianFilterVertical(DATA,sigma)

nData = size(DATA,1);
halfWindowSize = round(3*sigma);
windowSize = 2*halfWindowSize+1;
midpoint = halfWindowSize+1;
X = (1:windowSize)-midpoint;
Y = exp(-0.5*X.^2./sigma^2)/(sqrt(2*pi)*sigma);
C = conv2(DATA,Y');
indices2return = midpoint:(nData+midpoint);
F = conv2(DATA,Y','same');
end
