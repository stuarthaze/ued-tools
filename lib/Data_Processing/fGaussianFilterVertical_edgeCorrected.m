function F = fGaussianFilterVertical_edgeCorrected(DATA,sigma)

nData = size(DATA,1);
halfWindowSize = round(3*sigma);
windowSize = 2*halfWindowSize+1;
midpoint = halfWindowSize+1;
X = (1:windowSize)-midpoint;
Y = exp(-0.5*X.^2./sigma^2)/(sqrt(2*pi)*sigma);

DATA2 = zeros((nData+2*windowSize),1);
DATA2(1:windowSize,1) = ones(windowSize,1)*DATA(1);
dataIndicies = (windowSize+1):(windowSize+nData);
DATA2(dataIndicies) = DATA;
DATA2((windowSize+nData+1):(nData+2*windowSize)) = ones(windowSize,1)*DATA(nData);
% C = conv2(DATA,Y');
%C = conv2(DATA2,Y');
%F = C(indices2return,:);
F2 = conv2(DATA2,Y','same');
F = F2(dataIndicies);
end
