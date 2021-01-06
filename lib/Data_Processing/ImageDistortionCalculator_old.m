classdef ImageDistortionCalculator
    
    properties
        IMG_A
        IMG_Aorig
        IMG_B
        IMG_diff
        imgsizeA
        imgsizeB
        imageSizeDiff
        MASK
        distortionParams
        distortionParamsPreFit
        XX
        YY
        XX_AtoB
        YY_AtoB
        XXdist
        YYdist
        refineIntensityScale
        refineCenX
        refineCenY
        refineXshift
        refineYshift
        refineMag
        refineAspectRatio
        refineAspectTheta
        refineRotation
        refineRadialQuad
        refineRadialCubic
        refineBaseline
        refiningParamsLogical
        parVals
        numRefining
    
    end %Properties
    
    methods
        function obj = ImageDistortionCalculator(imgA)
            obj.IMG_A = imgA;
            obj.IMG_Aorig = imgA;
            obj.imgsizeA = size(imgA);
            [obj.XX, obj.YY] = meshgrid(1:size(imgA,2),1:size(imgA,1));
            obj = obj.refineNone();
            obj.MASK = ones(size(imgA));
        end
        function obj = refineNone(obj)
            obj.refineIntensityScale = 0;
            obj.refineCenX = 0;
            obj.refineCenY = 0;
            obj.refineXshift = 0;
            obj.refineYshift = 0;
            obj.refineMag = 0;
            obj.refineAspectRatio = 0;
            obj.refineAspectTheta = 0;
            obj.refineRotation = 0;
            obj.refineRadialQuad = 0;
            obj.refineRadialCubic = 0;
            obj.refineBaseline = 0;
        end
        function obj = loadIMG_B(obj,imgB)
            obj.IMG_B = imgB;
            obj.imgsizeB = size(imgB);
        end
        function obj = interpolateSizeAtoB(obj)
            [XXb,YYb] = meshgrid(1:obj.imgsizeB(2),1:obj.imgsizeB(1));
            XXab = XXb*obj.imgsizeA(2)/obj.imgsizeB(2);
            YYab = YYb*obj.imgsizeA(1)/obj.imgsizeB(1);
            obj.IMG_A = interp2(obj.XX,obj.YY,obj.IMG_A,XXab,YYab);
            obj.IMG_Aorig = obj.IMG_A;
            obj.imgsizeA = size(obj.IMG_A);
            obj.XX = XXb;
            obj.YY = YYb;
            obj.XX_AtoB = XXab;
            obj.YY_AtoB = YYab;
        end
        function obj = loadMask(obj,mask)
            obj.MASK = mask;
        end
        function obj = loadParams(obj,params)
            obj.distortionParams = params;
        end
        function obj = calcDistortion(obj)
            [obj.XXdist,obj.YYdist] = calcImageDistortion(obj.XX,obj.YY,obj.distortionParams);
        end
        function obj = calcIMG_A(obj)
            obj = obj.calcDistortion();
            img = calcDistortedIMG(obj.XX,obj.YY,obj.IMG_Aorig,obj.XXdist,obj.YYdist);
            obj.IMG_A = obj.distortionParams.intensityScale*img + obj.distortionParams.baseline;
        end
        function obj = updateRefiningParams(obj)
            rp(1) = obj.refineIntensityScale;
            rp(2) = obj.refineCenX;
            rp(3) = obj.refineCenY;
            rp(4) = obj.refineXshift;
            rp(5) = obj.refineYshift;
            rp(6) = obj.refineMag;
            rp(7) = obj.refineAspectRatio;
            rp(8) = obj.refineAspectTheta;
            rp(9) = obj.refineRotation;
            rp(10) = obj.refineRadialQuad;
            rp(11) = obj.refineRadialCubic;
            rp(12) = obj.refineBaseline;
            obj.refiningParamsLogical = logical(rp);
            obj.numRefining = sum(rp);
            if isempty(obj.parVals) && ~isempty(obj.distortionParams)
                obj.parVals = struct2parVals(obj.distortionParams);
            end
        end

        function PAR = refiningPars2allPars(obj,par)
            PAR = obj.parVals;
            PAR(obj.refiningParamsLogical) = par;
        end
        function IMGout = refiningPars2IMG(obj,pars)
            parsAll = obj.refiningPars2allPars(pars);
            parStruct = parVals2struct(parsAll);
            [X2,Y2] = calcImageDistortion(obj.XX,obj.YY,parStruct);
            IMGout = parsAll(1)*calcDistortedIMG(obj.XX,obj.YY,obj.IMG_Aorig,X2,Y2)+parsAll(12);
        end
        function obj = fitDistortionParams(obj)
            if ~((obj.imgsizeA(1)==obj.imgsizeB(1)) && (obj.imgsizeA(2)==obj.imgsizeB(2)))
                disp('Error: image size difference');
                return;
            end
            obj = obj.updateRefiningParams;
            p0 = obj.parVals(obj.refiningParamsLogical);
            fun2minimize = @(X) sum(sum((obj.MASK.*(obj.IMG_B - obj.refiningPars2IMG(X))).^2));
%             options = optimset('Algorithm','quasi-newton','Display','iter','MaxFunEvals',1e3);
            options = optimset('Display','iter','MaxFunEvals',1e3);
            pFit = fminunc(@(pFit)fun2minimize(pFit),p0,options);
            p0
            pFit
            obj.parVals(obj.refiningParamsLogical) = pFit;
            obj.distortionParamsPreFit = obj.distortionParams;
            obj.distortionParams = parVals2struct(obj.parVals);
            obj = obj.calcIMG_A;
            obj.IMG_diff = obj.MASK.*(obj.IMG_A - obj.IMG_B);
        end
        function [Xorig,Yorig] = XYdist2orig(obj,Xin,Yin)
            nx = length(Xin);
            ny = length(Yin);
            if nx ~= ny
                disp('x and y must be equal length');
                return
            end
            Xundist = interp2(obj.XXdist,Xin,Yin);
            Yundist = interp2(obj.YYdist,Xin,Yin);
            Xorig = interp2(obj.XX_AtoB,Xundist,Yundist);
            Yorig = interp2(obj.YY_AtoB,Xundist,Yundist);
        end
    end %methods
end %class
function p = struct2parVals(STRUCT)
    p(1) = STRUCT.intensityScale;
    p(2) = STRUCT.cen.x;
    p(3) = STRUCT.cen.y;
    p(4) = STRUCT.xshift;
    p(5) = STRUCT.yshift;
    p(6) = STRUCT.magnification;
    p(7) = STRUCT.aspectRatio;
    p(8) = STRUCT.aspectTheta;
    p(9) = STRUCT.rotation;
    p(10)= STRUCT.radialQuad;
    p(11)= STRUCT.radialCubic;
    p(12)= STRUCT.baseline;
end
function STRUCT = parVals2struct(p)
    STRUCT.intensityScale = p(1);
    STRUCT.cen.x = p(2);
    STRUCT.cen.y = p(3);
    STRUCT.xshift = p(4);
    STRUCT.yshift = p(5);
    STRUCT.magnification = p(6);
    STRUCT.aspectRatio = p(7);
    STRUCT.aspectTheta = p(8);
    STRUCT.rotation = p(9);
    STRUCT.radialQuad = p(10) ;
    STRUCT.radialCubic = p(11);
    STRUCT.baseline = p(12);
end
function IMG2 = calcDistortedIMG(XX1,YY1,IMG1,XX2,YY2)
    IMG2 = interp2(XX1,YY1,IMG1,XX2,YY2);
    IMG2(isnan(IMG2)) = 0;
end