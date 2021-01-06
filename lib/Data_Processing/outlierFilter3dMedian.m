function [IMGsOUT,outliers3d] = outlierFilter3dMedian(IMGs3d,threshold)
    [ny,nx,nz] = size(IMGs3d);
    if nz == 1
        [IMGsOUT,outliers3d] = outlierFilter2dMedian(IMGs3d,threshold);
    else
        IMGsOUT = IMGs3d;
        outliers3d = zeros([ny,nx,nz],'logical');
        for z = 1:nz
            [IMGsOUT(:,:,z),outliers3d(:,:,z)] = outlierFilter2dMedian(IMGs3d(:,:,z),threshold);
        end
    end
end