function out = displayIMGs_v2(IMGs,nRows,nCols,axisLimits,scaleMinMax,median2dFilter,gaussian2dFilter,timeDelays)
% axisLimits is a vector: [xmin, xmax, ymin, ymax]

nz = size(IMGs,3);
numberFigures =ceil(nz/(nRows*nCols));
imgsPerFig = nRows*nCols;

for fig = 1:numberFigures
    figure();
    for row = 1:nRows
        for col = 1:nCols
            iImg = (fig-1)*imgsPerFig + (row-1)*nCols + col;
            if iImg <= nz
                subplot(nRows,nCols,(row-1)*nCols + col);
                img = IMGs(:,:,iImg);
                if (nargin > 4) && median2dFilter
                    img = medfilt2(img,[median2dFilter,median2dFilter]);
                end
                if (nargin > 5) && gaussian2dFilter
                    img = imgaussfilt(img,gaussian2dFilter);
                end
                imagesc(img,scaleMinMax);
                if ~isempty(axisLimits)
                    axis(axisLimits);
                    nx = axisLimits(2)-axisLimits(1);
                    ny = axisLimits(4)-axisLimits(3);
                else
                    nx = size(img,2);
                    ny = size(img,1);
                end
                pbaspect([nx,ny,1])
                colormap('gray');
                if nargin > 6
                    title(num2str(timeDelays(iImg)));
                else
                    title(num2str(iImg));
                end
            else
                return
            end
        end
    end
end