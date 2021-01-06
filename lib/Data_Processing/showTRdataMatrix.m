function showTRdataMatrix(DATA,scalingQuantile)
    load('MyColormap_BWR','Colormap_BlueWhiteRed'); 
    ndat = size(DATA,1)*size(DATA,2);
    scale = quantile(reshape(abs(DATA),[ndat,1]),scalingQuantile);
    SCALE = [-scale,scale];
    figname = figure(); imagesc(DATA,SCALE);
    set(figname,'Colormap',Colormap_BlueWhiteRed);
    xlabel('Peak Index'); ylabel('Time Index');
    colorbar();
end