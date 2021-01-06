function [BG_2d, Radial_BG] = radialQuantile_BG_nestedEdit(IMG, centre, Qnt)
% IMG is a Matlab array
% centre is a structure with elements centre.x and centre.y
% Qnt is the quantile to use as the background estimate (0.5 = median value)

MASK = zeros(size(IMG));
Radial_quantile = fRadialQuantile(IMG,centre,MASK,Qnt);
R = 1:length(Radial_quantile);
RRadial = R.*Radial_quantile;
RRadial_NaNedit = replaceNaNsWithNearestValue(RRadial);
RRadial_smoothed = fGaussianFilterVertical(RRadial_NaNedit',1.5)';
Radial_BG = RRadial_smoothed./R;
IMG_Rmat = IMG_calcDists2Centre(IMG,centre);
IMG_Rmat(IMG_Rmat < 1) = 1;
BG_2d = Radial_BG(round(IMG_Rmat));

	function RadialQuantile = fRadialQuantile(IMG,CEN,MASK,Q)
	nx = size(IMG,2);
	ny = size(IMG,1);
	dmax(1) = sqrt(CEN.x^2+CEN.y^2);
	dmax(2) = sqrt((CEN.x-nx)^2+CEN.y^2);
	dmax(3) = sqrt(CEN.x^2+(ny-CEN.y)^2);
	dmax(4) = sqrt((CEN.x-nx)^2+(ny-CEN.y)^2);
	maxPix = round(max(dmax));
	intensities{maxPix} = [];
	for x = 1:nx
		for y = 1:ny
			if ~MASK(y,x)
				distSqrd = (x-CEN.x)^2 + (y-CEN.y)^2;
				r = round(sqrt(distSqrd));
				if distSqrd ~= 0
					n = length(intensities{r});
					intensities{r}(n+1) = IMG(y,x);
				end
			end
		end
	end
	nr = length(intensities)
	for RR = 1:nr
		RadialQuantile(RR) = quantile(intensities{RR},Q);
    end
    end

	function DATA2 = replaceNaNsWithNearestValue(DATA)
	nData = length(DATA);
	X = 1:nData;
	DATA2 = DATA;
	
	NANs = isnan(DATA);
	NANindicies = X(NANs);
	notNaNs = ~NANs;
	notNaNindicies = X(notNaNs);
	numNaNs = sum(NANs);

	if numNaNs
		for y = 1:numNaNs
			IndxDiff = abs(NANindicies(y)-notNaNindicies);
			[value, indx] = min(IndxDiff);
			DATA2(NANindicies(y)) = DATA(notNaNindicies(indx));
		end
	end
	end
	
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
	
	function R_IMG = IMG_calcDists2Centre(testIMG,centre)
	[ny,nx] = size(testIMG);
	Y = (1:ny)-centre.y;
	X = (1:nx)-centre.x;
	[XX,YY] = meshgrid(X,Y);
	R_IMG = sqrt(XX.^2 + YY.^2);
	end
	
end