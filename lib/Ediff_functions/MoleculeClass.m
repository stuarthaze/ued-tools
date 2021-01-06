classdef MoleculeClass
    
    properties
        XYZ
        XYZ_TR %Sequence of structures
        nAt
        atomicNumbers
        distList
        distList_TR
        nDists
        atomPairs
        covalentRadii
        bondedDists_Logical
        bondedPairs
        bondLengths
        bondLengths_TR
        nBonds
        numberingInCrystal
        bondedPairsCrystal

    end
    methods
        %Constructor
        function obj = MoleculeClass(atomicNumbers,XYZ)
            if nargin == 2
                obj.atomicNumbers = atomicNumbers;
                obj.XYZ = XYZ;
                obj.nAt = size(obj.XYZ,1);
                obj = obj.findBondedAtoms(1.3);
            end
        end
        
        function obj = readGeomFromCrystal(obj,CRYSTAL,atomIndices)           
            obj.atomicNumbers = CRYSTAL.atomicNumbers(atomIndices);
            obj.numberingInCrystal = atomIndices;
            obj.XYZ = CRYSTAL.XYZ(atomIndices,:);
            obj.nAt = size(obj.XYZ,1);
            obj = obj.findBondedAtoms(1.3);
        end
        
        function obj = findCovalentRadii(obj)
            obj.covalentRadii = readCovalentRadiiFromFile(obj.atomicNumbers);
        end
        
        function obj = createAtomPairList(obj)
            nAtms = obj.nAt;
            obj.nDists = nAtms*(nAtms-1)/2;
            obj.atomPairs = zeros([obj.nDists,2]);
            iDist = 0;
            for A = 1:(nAtms-1)
                for B = (A+1):nAtms
                    iDist = iDist+1;
                    obj.atomPairs(iDist,:) = [A,B];
                end
            end
        end
        
        function obj = calcDistList(obj)
            if isempty(obj.atomPairs)
                obj = obj.createAtomPairList();
            end
            obj.distList = calcDists(obj.atomPairs,obj.XYZ);
        end
        
        function obj = findBondedAtoms(obj,bondingThreshold)
            if isempty(obj.distList)
                obj = obj.calcDistList();
            end
            obj = obj.findCovalentRadii();
            sumCovalentRadii = sum(obj.covalentRadii(obj.atomPairs),2);
            obj.bondedDists_Logical = obj.distList <= bondingThreshold*sumCovalentRadii;
            obj.bondedPairs = obj.atomPairs(obj.bondedDists_Logical,:);
            obj.bondLengths = obj.distList(obj.bondedDists_Logical);
            obj.nBonds = sum(obj.bondedDists_Logical);
            if ~isempty(obj.numberingInCrystal)
                obj.bondedPairsCrystal = obj.numberingInCrystal(obj.bondedPairs);
            end
        end
        
        function obj = updateGeom(obj,XYZnew)
            obj.XYZ = XYZnew;
            obj = obj.calcDistList();
            obj.bondLengths = obj.distList(obj.bondedDists_Logical);
        end
        
        function bonds = xyz2bondlengths(obj,XYZnew)
            obj = obj.updateGeom(XYZnew);
            bonds = obj.bondLengths;
        end
        
%         function distsTR = calcBondLengthsTR(obj,XYZT,useCrystalNumbering,plotResults)
%             if nargin >= 3 && useCrystalNumbering
%                 atPairs = obj.numberingInCrystal(obj.bondedPairs);
%             else
%                 atPairs = obj.bondedPairs;
%             end
%             distsTR = calcDists(atPairs,XYZT);
%             if nargin == 4 && plotResults
%                 figure();
%                 plot(distsTR');
%                 for bond = 1:obj.nBonds
%                     legendtext{bond} = [num2str(atPairs(bond,1)),' - ',num2str(atPairs(bond,2))];
%                 end
%                 legend(legendtext);
%             end
%         end
%         
        function obj = readGeom_TR(obj,XYZT,fromCrystal)
            if nargin >= 3 && fromCrystal
                obj.XYZ_TR = XYZT(obj.numberingInCrystal,:,:);
            else
                obj.XYZ_TR = XYZT;
            end
        end

        function obj = calcBondLengths_TR(obj,XYZT,fromCrystal)
            obj = obj.readGeom_TR(XYZT,fromCrystal);
            obj.distList_TR = calcDists(obj.atomPairs,obj.XYZ_TR);
            obj.bondLengths_TR = obj.distList_TR(obj.bondedDists_Logical,:);
        end
        
        function plotBondLengths_TR(obj,useCrystalNumbering)
            if nargin > 1 && useCrystalNumbering
                atPairs = obj.bondedPairsCrystal;
            else
                atPairs = obj.bondedPairs;
            end
            plot(obj.bondLengths_TR');
            hold on
            for bond = 1:obj.nBonds
                refline(0,obj.bondLengths(bond));
                legendtext{bond} = [num2str(atPairs(bond,1)),' - ',num2str(atPairs(bond,2))];
            end
            legend(legendtext);
        end
        
    end %methods
end %classdef
%--------------------------------------------------------------
function dists = calcDists(atomPairs,XYZ)
    if ndims(XYZ) == 3
        dXYZ = XYZ(atomPairs(:,1),:,:) - XYZ(atomPairs(:,2),:,:);
        dists = sqrt(sum(dXYZ.^2,2));
        dists = reshape(dists,[size(dists,1),size(dists,3)]);
    else
        dXYZ = XYZ(atomPairs(:,1),:) - XYZ(atomPairs(:,2),:);
        dists = sqrt(sum(dXYZ.^2,2));
    end
end
