function M2 = expandMaskEdges(M1,npix)

Mdouble = double(M1);
K = createCircle(npix);
M2double = conv2(Mdouble,K,'same');
M2 = logical(M2double);

    function C = createCircle(radius)
        X = -radius:radius;
        [XX,YY] = meshgrid(X,X);
        RR = sqrt(XX.^2 + YY.^2);
        C = RR <= radius;
        C = double(C);
    end
end