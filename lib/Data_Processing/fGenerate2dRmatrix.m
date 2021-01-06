function Rmat = fGenerate2dRmatrix(cenX,cenY,MatSize);
[XX,YY] = meshgrid(1:MatSize(2),1:MatSize(1));
Rmat = sqrt((XX-cenX).^2 + (YY-cenY).^2);
% Rmat = zeros(MatSize);
% for Y = 1:MatSize(1)
%     for X = 1:MatSize(2)
%         Rmat(Y,X) = sqrt((X-cenX)^2 + (Y-cenY)^2);
%     end
% end
end