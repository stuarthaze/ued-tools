function [IMGsOUT,outliers3d] = outlierFilter3dGauss(IMGs3d,threshold)
    [ny,nx,nz] = size(IMGs3d);
    if nz == 1
        [IMGsOUT,outliers3d] = outlierFilter2dGauss(IMGs3d,threshold);
    else
        IMGsOUT = IMGs3d;
        outliers3d = zeros([ny,nx,nz],'logical');
        for z = 1:nz
            [IMGsOUT(:,:,z),outliers3d(:,:,z)] = outlierFilter2dGauss(IMGs3d(:,:,z),threshold);
        end
    end
end