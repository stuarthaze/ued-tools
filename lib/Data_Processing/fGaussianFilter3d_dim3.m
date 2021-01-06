function F = fGaussianFilter3d_dim3(DATA,sigma)
[nY,nX,nZ] = size(DATA);
F = zeros(size(DATA));
for y = 1:nY
    for x = 1:nX
        F(y,x,:) = fGaussian1d(squeeze(DATA(y,x,:)),sigma);
    end
end
    function F1d = fGaussian1d(D,sig)
        nData = length(D);
        halfWindowSize = round(3*sig);
        windowSize = 2*halfWindowSize+1;
        midpoint = halfWindowSize+1;
        X = (1:windowSize)-midpoint;
        Y = exp(-0.5*X.^2./sig^2)/(sqrt(2*pi)*sig);
        C = conv2(D,Y);
        F1d = conv2(D,Y,'same');
    end
end
