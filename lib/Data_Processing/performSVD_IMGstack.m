function results = performSVD_IMGstack(IMGs)
% (c) Stuart Hayes 23-01-2020
% Input: 3d array of nz images [ny,nx,nz]
% Output: 
%   results.imgs [ny,nx,nz]
%   results.svals [nz,1]
%   results.vmat [nz,nz]

[ny,nx,nz] = size(IMGs);

D = reshape(IMGs,[nx*ny,nz]);
[U,S,V] = svd(D);

results.imgs = reshape(U(:,1:nz),[ny,nx,nz]);
results.svals = diag(S);
results.vmat = V;

end
