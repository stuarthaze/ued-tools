classdef RadialDiffractionSimulator
    
    properties
        
        Iradial
        qvec
        dq_rad
        nx
        ny
        xcen
        ycen
        pixelSize
        distance
        wavelength
        astigmatism
        angleAstig
        magnification
        quadraticDist
        cubicDist
        intensityScale
        qIMG
        XX
        YY
        RR
        IMG
        parTable
        
    end %properties
    
    methods
        function obj = RadialDiffractionSimulator(NX,NY)
            obj.nx = NX;
            obj.ny = NY;
            [obj.XX, obj.YY] = meshgrid(1:NX,1:NY);
            obj.pixelSize = 24e-6;
            obj.distance = 0.5;
            obj.astigmatism = 1; % no distortion
            obj.angleAstig = 0;
            obj.magnification = 1;
            obj.quadraticDist = 0;
            obj.cubicDist = 0;
            obj.xcen = 0.5*NX;
            obj.ycen = 0.5*NY;
            obj.intensityScale = 1;
        end
        
        function obj = calcDistortion(obj)
            [X1,Y1] = meshgrid(1:obj.nx,1:obj.ny);
            dX = X1 - obj.xcen;
            dY = Y1 - obj.ycen;
            [dX2,dY2] = astigmaticDistortion(dX,dY,obj.astigmatism,obj.angleAstig);
            dRR = sqrt(dX2.^2 + dY2.^2);
            dX = dX2.*(1 + 1e-3*dRR*obj.quadraticDist + 1e-6*dRR.^2*obj.cubicDist)*obj.magnification;
            dY = dY2.*(1 + 1e-3*dRR*obj.quadraticDist + 1e-6*dRR.^2*obj.cubicDist)*obj.magnification;
            obj.XX = obj.xcen + dX;
            obj.YY = obj.ycen + dY;
            obj.RR = sqrt(dX.^2 + dY.^2);
        end
        
        function obj = calc_qIMG(obj)
            dX = obj.XX - obj.xcen;
            dY = obj.YY - obj.ycen;
            R = obj.pixelSize * sqrt(dX.^2 + dY.^2);
            THETA = atan(R./obj.distance);
            obj.qIMG = 4*pi*sin(THETA*0.5)/obj.wavelength;
        end
        
        function obj = calcIMG(obj)
            obj = obj.calcDistortion();
            obj = obj.calc_qIMG();
            obj.IMG = obj.intensityScale * interp1(obj.qvec,obj.Iradial,obj.qIMG);
            obj.IMG(isnan(obj.IMG)) = 0;
        end
        
        function obj = createParTable(obj)
            Names = {'xcen';'ycen';'astigmatism';'angleAstig';'magnification';'quadraticDist';'cubicDist';'intensityScale'};
            Values = [obj.xcen; obj.ycen; obj.astigmatism; obj.angleAstig; obj.magnification; obj.quadraticDist; obj.cubicDist; obj.intensityScale];
            Refine = logical(zeros([length(Values),1]));
            obj.parTable = table(Values,Refine);
            obj.parTable.Properties.RowNames = Names;
        end
        
        function obj = updateParTable(obj)
            obj.parTable.Values = [obj.xcen; obj.ycen; obj.astigmatism; obj.angleAstig; obj.magnification; obj.quadraticDist; obj.cubicDist; obj.intensityScale];
        end
        
        function obj = setParamsFromTable(obj)
            obj.xcen = obj.parTable.Values(1);
            obj.ycen = obj.parTable.Values(2);
            obj.astigmatism = obj.parTable.Values(3);
            obj.angleAstig = obj.parTable.Values(4);
            obj.magnification = obj.parTable.Values(5);
            obj.quadraticDist = obj.parTable.Values(6);
            obj.cubicDist = obj.parTable.Values(7);
            obj.intensityScale = obj.parTable.Values(8);
        end
        
        function vals = readRefiningParams(obj)
            vals = obj.parTable.Values(obj.parTable.Refine);
        end
        
        function obj = parVec2table(obj,X)
            obj.parTable.Values(obj.parTable.Refine) = X;
            obj = obj.setParamsFromTable();
        end
        
        function img = pars2img(obj,X)
            obj = obj.parVec2table(X);
            obj = obj.calcIMG();
            img = obj.IMG;
        end
                
    end %methods

end %classdef
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