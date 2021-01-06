function out = viewCircularROIs(IMG, xROIs, yROIs, dROI, textSize, intensityScale)

nROIs = length(xROIs);
figure();
imgScale = intensityScale*max(max(IMG)); 
imagesc(-IMG, imgScale*[-1,0]);
colormap('gray');
pbaspect([1,1,1]);
hold on;
rROI = dROI*0.5;
th = 0:0.1:2*pi;
x_circle = rROI*cos(th);
y_circle = rROI*sin(th);
for i = 1:nROIs
    x = xROIs(i) + x_circle;
    y = yROIs(i) + y_circle;
    plot(x,y,'Color','g')
    text(xROIs(i),yROIs(i),num2str(i),'FontSize',textSize,'Color',[0.4 0.4 1]);
end
hold off;