function [XXout,YYout] = calcImageDistortion(XXin,YYin,params)

    xx = XXin - params.cen.x;
    yy = YYin - params.cen.y;
    if params.aspectRatio ~= 1
        [xx,yy] = astigmaticDistortion(xx,yy,params.aspectRatio,params.aspectTheta);
    end
    if params.rotation
        [xx,yy] = rotateCoords(xx,yy,params.rotation);
    end
    if params.magnification ~= 1
        xx = xx/params.magnification;
        yy = yy/params.magnification;
    end
    if params.radialQuad || params.radialCubic
        [xx,yy] = radialDistortion(xx,yy,params.radialQuad,params.radialCubic);
    end
    if params.xshift
        xx = xx - params.xshift;
    end
    if params.yshift
        yy = yy - params.yshift;
    end
    
    XXout = xx + params.cen.x;
    YYout = yy + params.cen.y;

%----------------------------------------------
    function [XX2,YY2] = radialDistortion(XX1,YY1,quad,cubic)
    R2rd = XX1.^2 + YY1.^2;
    RR = sqrt(R2rd);
    XX2 = XX1.*(1 + quad*RR + cubic*R2rd);
    YY2 = YY1.*(1 + quad*RR + cubic*R2rd);
    end

    function [XX2,YY2] = rotateCoords(XX1,YY1,theta)
        XX2 = XX1*cos(theta) - YY1*sin(theta);
        YY2 = XX1*sin(theta) + YY1*cos(theta);
    end

    function [XX2,YY2] = astigmaticDistortion(XX1,YY1,ratio,theta)
        [XX2, YY2] = rotateCoords(XX1,YY1,theta);
        XX3 = XX2*ratio;
        YY3 = YY2/ratio;
        [XX2,YY2] = rotateCoords(XX3,YY3,-theta);
    end

end