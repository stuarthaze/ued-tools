classdef GEDcalculator
    % (c) SAH 2020
    
    properties
        atomicNumbers
        xyz
        nat
        nmol
        smin
        smax
        ds
        ns
        svals
        coherenceL
        Fatomic
        Iatomic
        Imol
        Itotal
        sMs
    end
    
    methods
        %Constructor
        function obj = GEDcalculator(atomicZ,XYZ)
            obj.atomicNumbers = atomicZ;
            obj.nat = length(atomicZ);
            obj = obj.loadxyz(XYZ);
        end
        
        function obj = loadxyz(obj,XYZ)
            if ndims(XYZ) == 3
                obj.nmol = size(XYZ,3);
            else
                obj.nmol = 1;
            end
            obj.xyz = XYZ;
        end
        
        function obj = setSvals(obj,SMIN,SMAX,DS)
            obj.smin = SMIN;
            obj.smax = SMAX;
            obj.ds = DS;
            obj.svals = SMIN:DS:SMAX;
            obj.ns = length(obj.svals);
        end
        
        function obj = generateFatomic(obj)
            KrkTbl = load('KirklandTable.mat','-ascii');
            obj.Fatomic = zeros([obj.nat,obj.ns]);
            for a = 1:obj.nat
                for b = 1:obj.ns
                    s = obj.svals(b);
                    obj.Fatomic(a,b) = KirklandTable2f(obj.atomicNumbers(a),s,KrkTbl);
                end
            end
        end
        
        %---------------------------------------
        function obj = calculateDiffraction(obj,coherenceLength)
            if nargin == 2
                obj.coherenceL = coherenceLength;
            end
            if isempty(obj.Fatomic)
                obj = obj.generateFatomic();
            end
            obj.Iatomic = sum(obj.Fatomic.^2,1);

            if obj.nmol == 1
                obj.Imol = calcImol_v3(obj.svals, obj.Fatomic, obj.xyz, obj.coherenceL);
            elseif obj.nmol > 1
                Imoln = zeros(obj.nmol,obj.ns);
                for mol = 1:obj.nmol
                    Imoln(mol,:) = calcImol_v3(obj.svals, obj.Fatomic, obj.xyz(:,:,mol),obj.coherenceL);
                end
                obj.Imol = mean(Imoln,1);
            end
            % Included in calcImol_v3
%             % Correct for division by zero at s=0
%             if obj.svals(1) == 0
%                 I0 = 0;
%                 for a = 1:(obj.nat-1)
%                     for b = (a+1):obj.nat
%                         I0 = I0 + obj.Fatomic(a,1)*obj.Fatomic(b,1);
%                     end
%                 end
%                 obj.Imol(1) = 2*I0;
% %                 obj.Imol(1) = 0;
%             end
            obj.Itotal = obj.Iatomic + obj.Imol;
            obj.sMs = obj.svals.*obj.Imol./obj.Iatomic;
        end
        
    end %methods
end %class


