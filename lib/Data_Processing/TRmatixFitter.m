classdef TRdataFitter

    properties
        % 2d arrays [ndat x nt]
        Mdata
        Mcalc
        Mdiff
        Tvec
        ndat
        nt
        nexp
        % Parameters
        refineBase
        Base
        refineSteps % 1 or 0
        AmpStep     %size=ndat
        refineExp   %size=nexp
        AmpExp      %size=[ndat x nexp]
        TauExp      %size=[ndat x nexp]
        refineT0    
        t0
        refineConv
        tConv
        
        parIndxBase
        parIndxStep
        parIndxAmpExp % ndat x nexp
        parIndxTau
        parIndxT0
        parIndxConv
        nPar
        nParRef
        refiningParsLogical
        parVals
        
        fittingOptions
        
    end %properties
    
    methods
        %Constructor
        function obj = TRmatixFitter(DATA,xvec)
            obj.Mdata = DATA;
            [obj.ndat,obj.nt] = size(DATA);
            obj.Tvec = xvec;
        end
        
        function obj = initializeParameters(obj)
            obj.Base = zeros([obj.ndat,1]);
            obj.refineBase = FALSE;
            obj.AmpStep = zeros([obj.ndat,1]);
            obj.refineSteps = FALSE;
            obj.refineExp = logical(zeros([1,obj.nexp]));
            obj.AmpExp = zeros([obj.ndat,obj.exp]);
            obj.TauExp = ones([1,obj.exp]);
        end
        
        function obj = calculateParIndices(obj)
            indx = 0;
            obj.parIndxBase = 1:obj.ndat;
            obj.parIndxStep = indx + 1:obj.ndat;
            for iExp = 1:obj.nexp
                if obj.refineExp(iExp)
                    obj.parIndxAmpExp(:,iExp) = indx + 1:obj.ndat;
                    indx = indx + obj.ndat;
                    obj.parIndxTau(iExp) = indx + 1;
                    indx = indx + 1;
                end
            end
            obj.parIndxT0 = indx + 1;
            indx = indx + 1;
            obj.parIndxConv = indx + 1;
            indx = indx + 1;
            obj.nPar = indx;
        end
        %---------------------------------------
        function obj = setRefiningParams(obj)
            obj.refiningParsLogical = logical(zeros([1,obj.nPar]));
            if obj.refineBase
                obj.refiningParsLogical(obj.parIndxBase) = 1;
            else
                obj.refiningParsLogical(obj.parIndxBase) = 0;
            end
            if obj.refineSteps
                obj.refiningParsLogical(obj.parIndxStep) = 1;
            else
                obj.refiningParsLogical(obj.parIndxStep) = 0;
            end
            for iExp = 1:obj.nexp
                if obj.refineExp(iExp)
                    obj.refiningParsLogical(obj.parIndxAmpExp(:,iExp)) = 1;
                    obj.refiningParsLogical(obj.parIndxTau(iExp)) = 1;
                else
                    obj.refiningParsLogical(obj.parIndxAmpExp(:,iExp)) = 0;
                    obj.refiningParsLogical(obj.parIndxTau(iExp)) = 0;
                end
            end
            if obj.refineT0
                obj.refiningParsLogical(obj.parIndxT0) = 1;
            else
                obj.refiningParsLogical(obj.parIndxT0) = 0;
            end
            if obj.refineConv
                obj.refiningParsLogical(obj.parIndxConv) = 1;
            else
                obj.refiningParsLogical(obj.parIndxConv) = 0;
            end
            obj.nParRef = sum(obj.refiningParsLogical);
        end
        
        function obj = generateParList(obj)
            par = zeros([1,obj.nPar]);
            par(obj.parIndxBase) = obj.Base;
            par(obj.parIndxStep) = obj.AmpStep;
            for iExp = 1:obj.nexp
                if obj.refineExp(iExp)
                    par(obj.parIndxAmpExp(:,iExp)) = obj.AmpExp(:,iExp);
                    par(obj.parIndxTau(iExp)) = obj.TauExp(iExp);
                end
            end
            par(obj.parIndxT0) = obj.t0;
            par(obj.parIndxConv) = obj.tConv;
            obj.parVals = par;
        end
        
        function obj = updateObjectPars(obj)
            obj.Base = obj.parVals(obj.parIndxBase);
            obj.AmpStep = obj.parVals(obj.parIndxStep);
            for iExp = 1:obj.nexp
                if obj.refineExp(iExp)
                    obj.AmpExp(:,iExp) = obj.parVals(obj.parIndxAmpExp(:,iExp));
                    obj.TauExp(iExp) = obj.parVals(obj.parIndxTau(iExp));
                end
            end
            obj.t0 = obj.parVals(obj.parIndxT0);
            obj.tConv = obj.parVals(obj.parIndxConv);
        end
        
        
        %------------------------------------
        function Mout = params2mcalc(obj,PARsIn)
            par = obj.parVals;
            par(obj.refiningParsLogical) = PARsIn;
            %set variable values
            T0 = par(obj.parIndxT0);
            base = par(obj.parIndxBase);
            Astep = par(obj.parIndxStep);
            for iExp = 1:obj.nexp
                if obj.refineExp(iExp)
                    Aexp(:,iExp) = par(obj.parIndxAmpExp(:,iExp));
                    Texp(iExp) = par(obj.parIndxTau(iExp));
                end
            end
            tCon = par(obj.parIndxConv);
                
            dT = obj.Tvec - T0;
            Tneg = dt <= 0;
            
            stepT = 0.5*(1+erf(dT/(tCon*sqrt(2))));
            expT = zeros([obj.nexp,obj.nt]);
            
            for ii = 1:obj.nexp
                unitExp = 1 - exp(-dT/Texp(ii));
                unitExp(Tneg) = 0;
                expT(ii,:) = imgaussfilt(unitExp,tCon);
            end
            
            Mout = zeros([obj.ndat,obj.nt]);
            for dat = 1:obj.ndat
                Mline = base(dat)+ Astep(dat)*stepT;
                for iExp = 1:obj.nexp
                    Mline = Mline + Aexp(dat,iExp)*expT(iExp,:);
                end
                Mout(dat,:) = Mline;
            end
        end %params2mcalc
        %------------------------------------
        % Fitting
        
        
    end %methods
    
end %classdef
