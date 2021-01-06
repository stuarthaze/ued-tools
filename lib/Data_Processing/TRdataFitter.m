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
        parValsPrevious
                
    end %properties
    
    methods
        %Constructor
        function obj = TRdataFitter(DATA,xvec)
            obj.Mdata = DATA;
            [obj.ndat,obj.nt] = size(DATA);
            obj.Tvec = xvec;
        end
        
        function obj = initializeParameters(obj)
            obj.Base = zeros([obj.ndat,1]);
            obj.refineBase = 0;
            obj.AmpStep = zeros([obj.ndat,1]);
            obj.refineSteps = 0;
            obj.refineExp = zeros([1,obj.nexp]);
            obj.AmpExp = zeros([obj.ndat,obj.nexp]);
            obj.TauExp = ones([1,obj.nexp]);
            obj = obj.calculateParIndices();
        end
        
        function obj = calculateParIndices(obj)
            indx = 0;
            obj.parIndxBase = 1:obj.ndat;
            indx = obj.ndat;
            obj.parIndxStep = indx + (1:obj.ndat);
            indx = indx + obj.ndat;
            for iExp = 1:obj.nexp
                obj.parIndxAmpExp(1:obj.ndat,iExp) = indx + (1:obj.ndat);
                indx = indx + obj.ndat;
                obj.parIndxTau(iExp) = indx + 1;
                indx = indx + 1;
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
            if length(obj.refineExp) < obj.nexp
                disp('Error in size(refineExp): setting size = nexp');
                obj.refineExp = obj.refineExp*ones([1,obj.nexp]);
            end
            if length(obj.TauExp) < obj.nexp
                disp('Error: Inconsistent number of time constants');
            end
            for iExp = 1:obj.nexp
                if obj.refineExp(iExp)
                    obj.refiningParsLogical(obj.parIndxAmpExp(1:obj.ndat,iExp)) = 1;
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
                par(obj.parIndxAmpExp(1:obj.ndat,iExp)) = obj.AmpExp(1:obj.ndat,iExp);
                par(obj.parIndxTau(iExp)) = obj.TauExp(iExp);
            end
            par(obj.parIndxT0) = obj.t0;
            par(obj.parIndxConv) = obj.tConv;
            obj.parVals = par;
        end
        
        function obj = updateObjectPars(obj)
            obj.Base = obj.parVals(obj.parIndxBase);
            obj.AmpStep = obj.parVals(obj.parIndxStep);
            for iExp = 1:obj.nexp
                obj.AmpExp(:,iExp) = obj.parVals(obj.parIndxAmpExp(:,iExp));
                obj.TauExp(iExp) = obj.parVals(obj.parIndxTau(iExp));
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
            Tneg = dT <= 0;
            dtDat = (obj.Tvec(end)-obj.Tvec(1))/(length(obj.Tvec)-1);
            
            stepT = 0.5*(1+erf(dT/(tCon*sqrt(2))));
            expT = zeros([obj.nexp,obj.nt]);
            
            for ii = 1:obj.nexp
                unitExp = 1 - exp(-dT/Texp(ii));
                unitExp(Tneg) = 0;
                expT(ii,:) = imgaussfilt(unitExp,tCon/dtDat);
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
        function obj = fitData(obj,fitOptions)
            p0 = obj.parVals(obj.refiningParsLogical);
%             
%             if obj.ndat == 1
%                 fun2fit = @(X) sum((obj.Mdata-obj.params2mcalc(X)).^2);
%             else
%                 fun2fit = @(X) sum(sum((obj.Mdata-obj.params2mcalc(X)).^2));
%             end
%             
%             pFit = fminunc(@(pFit)fun2fit(pFit),p0,fitOptions);

            fun2fit = @(X) obj.Mdata-obj.params2mcalc(X);
            lb = [];
            ub = [];
            pFit = lsqnonlin(@(pFit)fun2fit(pFit),p0,lb,ub,fitOptions);
            
            obj.parValsPrevious = obj.parVals;
            obj.parVals(obj.refiningParsLogical) = pFit;
            obj = obj.updateObjectPars();
            obj.Mcalc = obj.params2mcalc(pFit);
            obj.Mdiff = obj.Mdata - obj.Mcalc;

        end %fitting
        
        function showFit(obj)
            figure();
            if obj.ndat == 1
                plot(obj.Tvec,obj.Mdata,obj.Tvec,obj.Mcalc);
                legend('Data','Fit');
                xlabel('Time');
            else
                subplot(1,3,1); imagesc(obj.Mdata); colorbar(); title('Data');
                subplot(1,3,2); imagesc(obj.Mcalc); colorbar(); title('Fit');
                subplot(1,3,3); imagesc(obj.Mdiff); colorbar(); title('Difference');
                
            end
        end
        
    end %methods
    
end %classdef
