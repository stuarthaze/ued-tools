function OUT = saveAnnotatedMovieFrames(IMGstack,SCALE_MINMAX,timeVector,dirSave)

[ny,nx,nt] = size(IMGstack);
if ~isdir(dirSave)
    mkdir(dirSave);
end

figname = figure();
for T = 1:nt
    imagesc(IMGstack(:,:,T),SCALE_MINMAX);
    pbaspect([nx,ny,1]);
	colormap('gray');
    timetext = ['image ',num2str(T),'  ', num2str(round(timeVector(T))),' fs'];
    text(10,10,timetext);
    f = getframe(gcf);
    newimg = f.cdata;
    outputFilename = [dirSave,'\',timetext,'.bmp'];
    imwrite(newimg,outputFilename, 'bmp');
end
close(gcf);
end