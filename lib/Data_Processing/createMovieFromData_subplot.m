function createMovieFromData_subplot(IMGs,IMGs2,tDelays,SCALE,AxisLimits,frameRate,fileSaveName,txtPosition,title1,title2)


[ny,nx,nFrames] = size(IMGs);
if isempty(AxisLimits)
    AxisLimits = [1,nx,1,ny];
end
nxcrop = AxisLimits(2) - AxisLimits(1) + 1;
nycrop = AxisLimits(4) - AxisLimits(3) + 1;
xmin = AxisLimits(1);
ymin = AxisLimits(3);


% figure();
% ax = gca;
% ax.NextPlot = 'replaceChildren';
% mov(1:nFrames) = struct('cdata', [],'colormap', []);

h = figure();


scale2 = -mean(min(min(IMGs2,[],2),[],1));
SCALE2 = scale2*[-1,1];

for T = 1:nFrames
    subplot(1,2,1); imagesc(IMGs(:,:,T),SCALE);
    colormap('gray');
    title(title1);
%     if ~isempty(AxisLimits)
    axis(AxisLimits);
    pbaspect([nxcrop,nycrop,1]);
%     else
%         pbaspect([nx,ny,1]);
%     end
    set(gca, 'Ydir', 'reverse');
    if ~isempty(tDelays)
        txt = [num2str(tDelays(T)),' fs'];
        if ~isempty(txtPosition)
            dx = txtPosition(2);
            dy = txtPosition(1);
        else
            dx = 0;
            dy = 0;
        end
        t = text(xmin+dx,ymin+dy,txt);
        t.Color = 'white';
        t.FontSize = 8;
    end
%     title();
    subplot(1,2,2); imagesc(IMGs2(:,:,T),SCALE2);
    pbaspect([1,1,1]); title(title2);
    drawnow;
    mov(T) = getframe(gcf);
end


myVideo = VideoWriter([fileSaveName,'.avi']);
myVideo.FrameRate = frameRate;
myVideo.Quality = 100;
open(myVideo);
writeVideo(myVideo, mov);
close(myVideo);
end