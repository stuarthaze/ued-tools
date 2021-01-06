function createMovieFromData(IMGs,tDelays,SCALE,AxisLimits,frameRate,fileSaveName,txtPosition)


[ny,nx,nFrames] = size(IMGs);
if isempty(AxisLimits)
    AxisLimits = [1,nx,1,ny];
end
nxcrop = AxisLimits(2) - AxisLimits(1) + 1;
nycrop = AxisLimits(4) - AxisLimits(3) + 1;
xmin = AxisLimits(1);
ymin = AxisLimits(3);


figure();
ax = gca;
ax.NextPlot = 'replaceChildren';
mov(1:nFrames) = struct('cdata', [],'colormap', []);

for T = 1:nFrames
    imagesc(ax,IMGs(:,:,T),SCALE);
    colormap('gray');
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
    end
%     title();
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