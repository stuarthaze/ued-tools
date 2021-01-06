function viewImageWithCircles(IMG, intensityScaling, invertGrayscale, xCen, yCen, radii, textSize)
% Plots image with circles
% Set textSize = 0 to view without labels
% intensityScaling is a vector [min, max]
% invertGrayscale is a logical -> white is 0, black is maximum
% It would be good to make this more general in the future
% the number of circles is determined by the maximum of the number of radii or centres 
nx = length(xCen);
ny = length(yCen);
nr = length(radii);
if nx ~= ny
    disp('error in centres');
    return
end
if nx == nr
    xCenVec = xCen;
    yCenVec = yCen;
    rVec = radii;
elseif nr == 1
    xCenVec = xCen;
    yCenVec = yCen;
    rVec = radii*ones([1,nx]);
elseif nx == 1
    xCenVec = xCen*ones([1,nr]);
    yCenVec = yCen*ones([1,nr]);
    rVec = radii;
end
nCircles = length(xCenVec);

figure();
if invertGrayscale
    intensityScaling = -(fliplr(intensityScaling));
    IMG = -IMG;
    textColor = [0 0 1];
else
    textColor = [0 1 1];
end
imagesc(IMG, intensityScaling);
colormap('gray');
pbaspect([1,1,1]);
hold on;
th = (0:0.01:1)*2*pi;
x_circle = cos(th);
y_circle = sin(th);
for ii = 1:nCircles
    x = xCenVec(ii) + rVec(ii)*x_circle;
    y = yCenVec(ii) + rVec(ii)*y_circle;
    plot(x,y,'Color','g')
    if textSize
        text(xCenVec(ii)+rVec(ii),yCenVec(ii),num2str(ii),'FontSize',textSize,'Color',textColor);
    end
end
hold off;
end