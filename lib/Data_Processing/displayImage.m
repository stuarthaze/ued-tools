function displayImage(img,displayOptions)
% Fields for displayOptions:
%  gamma
%  scale
%  colormap
%  invert
% SAH 2020

    if (nargin == 2) && isfield(displayOptions,'gamma')
        img(img < 0) = 0;
        img = img.^displayOptions.gamma;
    end

    if (nargin == 2) && isfield(displayOptions,'scale')
        scale = displayOptions.scale;
        if length(scale) == 1
            scalemin = 0;
            scalemax = scale;
        elseif length(scale) == 2
            scalemin = scale(1);
            scalemax = scale(2);
        end
    else
        scalemin = min(min(img));
        scalemax = max(max(img));
    end

    % Display
    imagesc(img,[scalemin,scalemax]); 
    pbaspect([size(img,2),size(img,1),1]);


    if (nargin == 2) && isfield(displayOptions,'colormap')
        cmap = colormap(displayOptions.colormap);
    else
        cmap = colormap('gray');
    end

    if (nargin == 2) && isfield(displayOptions,'invert')
        colormap(flipud(cmap));
    end

end