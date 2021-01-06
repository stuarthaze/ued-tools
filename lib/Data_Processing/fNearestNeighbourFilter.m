function IMG2 = fNearestNeighbourFilter(IMG,nStdDevs )
% Function to set bad pixels to the averge of neighbours
% nStdDevs sets the number of sigma cut-off
imgSize = size(IMG);
IMG2(:,:) = IMG(:,:);
nx = imgSize(2);
ny = imgSize(1);
for x = 2:(nx-1)
    for y = 2:(ny-1)
        neighb = IMG([y-1:y+1],[x-1:x+1]);
        vals1d = [neighb(1:3),neighb(4),neighb(6),neighb(7:9)];
        av = mean(vals1d);
        err = std(vals1d);
        if (IMG(y,x) > av+nStdDevs*err) || (IMG(y,x) < av-nStdDevs*err)
            IMG2(y,x) = av;
        end
    end
end
end