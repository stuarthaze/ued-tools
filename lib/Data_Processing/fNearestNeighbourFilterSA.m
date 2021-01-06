function IMG2 = fNearestNeighbourFilterSA(IMG,nStdDevs,SAx,SAy)
% Function to set bad pixels to the averge of neighbours, selected area
% nStdDevs sets the number of sigma cut-off
% SAx, SAy - Areas to filter
imgSize = size(IMG);
IMG2(:,:) = IMG(:,:);
nx = imgSize(2);
ny = imgSize(1);
for x = SAx(1):SAx(2)
    for y = SAy(1):SAy(2)
        neighb = IMG([y-1:y+1],[x-1:x+1]);
        vals1d = [neighb(1:3),neighb(4),neighb(6),neighb(7:9)];
        av = mean(vals1d);
        err = std(vals1d);
        if (IMG(y,x) > av+nStdDevs*err) | (IMG(y,x) < av-nStdDevs*err)
            IMG2(y,x) = av;
        end
    end
end
end