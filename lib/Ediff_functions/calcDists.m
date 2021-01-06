function dists = calcDists(atomPairs,XYZ)
% Calculates distances from a list of atom pairs and cartesian coords
% Also handles a series of coordinates
    if ndims(XYZ) == 3
        dXYZ = XYZ(atomPairs(:,1),:,:) - XYZ(atomPairs(:,2),:,:);
        dists = sqrt(sum(dXYZ.^2,2));
        dists = reshape(dists,[size(dists,1),size(dists,3)]);
    else
        dXYZ = XYZ(atomPairs(:,1),:) - XYZ(atomPairs(:,2),:);
        dists = sqrt(sum(dXYZ.^2,2));
    end
end