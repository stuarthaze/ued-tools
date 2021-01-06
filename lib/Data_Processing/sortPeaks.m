function [vals_sorted, locs_sorted] = sortPeaks(PEAKS)

[ny,nx] = size(PEAKS);
ndata = ny*nx;
[vals_sorted, locs_sorted] = sort(reshape(PEAKS,1,ndata),'descend');
end