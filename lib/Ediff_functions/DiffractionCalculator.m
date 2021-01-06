classdef DiffractionCalculator
    
    properties
        crystal
        experiment
        model
        %Calculated internally
        nHKL
        HKL_list
        qHKL_list
        modQ_hkl
        Fatomic
        FatomicT
        xExcited
        UVW_GS
        UVW_ES
        Uiso_GS
        Uiso_ES
        Fhkl_GS
        Fhkl_ES
        Ihkl_GS
        Ihkl_ES
        Ihkl_diff
        
        nPix
        k0
        modk0
        qIMG_list
        maxqIMG
        modqIMG_list
        
        WTS %Matrix to weight reflections: [nPix x nHKL]
        IMG_GS
        IMG_diff
        
    end %properties
    %----------------------------------------------------------------------
    methods
        % Constructor
        function obj = DiffractionCalculator(CRYSTAL,EXPERIMENT,MODEL)
            obj.crystal = CRYSTAL;
            obj.experiment = EXPERIMENT;
            obj.model = MODEL;
            
            %--------------------------------------------------------
            %            Nested functions for constructor
%             
%             function Tmat = transformMatFrom_k0(k0)
%                 V3 = k0/norm(k0); % e-beam in z-direction
%                 Vy = cross(V3,[1,0,0]);
%                 V2 = Vy/norm(Vy);
%                 V1 = cross(V2,V3);
%                 Tmat = [V1;V2;V3];
%             end
            
            function [Q,Qmod] = calcQimg()
                phi = EXPERIMENT.imgRotation;
                Q = zeros([imgsize.y,imgsize.x,3]);
                [XX,YY] = meshgrid(1:imgsize.x,1:imgsize.y);
                dXX = XX-cen.x;
                dYY = YY-cen.y;
                RR = sqrt(dXX.^2 + dYY.^2);
                AA = atan(pixelscale*RR);
                Qmod = 2*modk0*sin(AA*0.5);
                Q(:,:,3) = -Qmod.*sin(AA);
                qRR = Qmod.*cos(AA);
                qXX = qRR.*dXX./RR;
                qYY = qRR.*dYY./RR;
                Q(:,:,1) = qXX*cos(phi)-qYY*sin(phi);
                Q(:,:,2) = qXX*sin(phi)+qYY*cos(phi);
            end
            %------------------
            % Calculate geometry in terms of q
            
            imgsize = EXPERIMENT.imgsize;
            cen = EXPERIMENT.cen;
            pixelscale = EXPERIMENT.pixelSize/(EXPERIMENT.detectorDist*EXPERIMENT.magnification); % angle between pixels at origin
            modk0 = 1/EXPERIMENT.wavelength;
            obj.experiment.dqPix = modk0*pixelscale;
            beamSizeSigma = EXPERIMENT.beamSizeFWHM/2.35;
            beamDivSigma = beamSizeSigma/(EXPERIMENT.detectorDist*EXPERIMENT.magnification);
            obj.experiment.beamDivSigma = beamDivSigma;
            obj.experiment.spotSizePixSigma = beamSizeSigma/EXPERIMENT.pixelSize;

            beamDir_xyz = EXPERIMENT.beamDir_uvw*CRYSTAL.T_FracCart;
            k0_xyz_norm = beamDir_xyz/norm(beamDir_xyz);
            k0 = k0_xyz_norm*modk0;
            k0_normT = k0_xyz_norm';
            [qIMG, modQ] = calcQimg();

            % Calculate hkl values to use
            nPix = imgsize.x*imgsize.y;
            sIMG_L = reshape(qIMG,[nPix,3]); % s-vectors for img
            T_SQ = transformMatFrom_k0(k0);
            qIMG_L = sIMG_L*T_SQ;
            modqIMG_L = reshape(modQ,[nPix,1]);

            HKL_IMG_L = round(qIMG_L*CRYSTAL.T_xyz2hkl);
            HKL_minvals = min(HKL_IMG_L);
            HKL_maxvals = max(HKL_IMG_L);
            HKL_maxList = createHKL_list(HKL_minvals,HKL_maxvals);
            qHKL_maxList = HKL_maxList*CRYSTAL.T_hkl2xyz;
            qHKLs_dotk0 = qHKL_maxList*k0_normT;
            qIMG_dotk0 = qIMG_L*k0_normT;
            nSigmaCutoff = 3;
            maxqIMG = max(modqIMG_L);
            maxSigmaParallel = sqrt((MODEL.thetaParallel*maxqIMG)^2 + (beamDivSigma*maxqIMG)^2 + MODEL.sigmaEwald^2);
            qCutoff = nSigmaCutoff*maxSigmaParallel
            qdotk0_min = min(qIMG_dotk0)-qCutoff
            qdotk0_max = max(qIMG_dotk0)+qCutoff
            modqHKL_maxList = sqrt(sum(qHKL_maxList.^2,2));
            indices2keep = (qHKLs_dotk0 >= qdotk0_min) & (qHKLs_dotk0 <= qdotk0_max) & (modqHKL_maxList <= maxqIMG);
            HKL_list = HKL_maxList(indices2keep,:);
            qHKL_list = qHKL_maxList(indices2keep,:);
            modQ_hkl = sqrt(sum(qHKL_list.^2,2));
            
            % Calculate Fatomic for HKL_list
            if strcmp(MODEL.fAtomicType,'Kirkland')
                Fatomic = calculateAtomicF_Kirkland(CRYSTAL,HKL_list);
            elseif strcmp(MODEL.fAtomicType,'EDICO')
                modS_hkl = modQ_hkl*2*pi;
                [F_mod,F_phi] = interp_f_from_tables_EDICO(CRYSTAL.atomicNumbers,modS_hkl,EXPERIMENT.voltage_kV);
                if MODEL.fAtomicComplex
                    Fatomic = (F_mod.*exp(1i*F_phi))';
                else
                    Fatomic = F_mod';
                end
            elseif ~exist('MODEL.fAtomicType','var')
                disp('Need to specify atomic scattering data source');
            end 

            %Set object properties required later
            obj.HKL_list = HKL_list;
            obj.nHKL = size(HKL_list,1);
            obj.qHKL_list = qHKL_list;
            obj.modQ_hkl = modQ_hkl;
            obj.k0 = k0;
            obj.modk0 = modk0;
            obj.nPix = nPix;
            obj.qIMG_list = qIMG_L;
            obj.maxqIMG = maxqIMG;
            obj.modqIMG_list = modqIMG_L;
            obj.Fatomic = Fatomic;
            obj.FatomicT = Fatomic';

        end %constructor
        
        %-----------------------------------------------------------

        function obj = calculateWTS(obj)
            obj.WTS = zeros([obj.nPix,obj.nHKL]);
            k0_mesh = meshgrid(obj.k0, 1:obj.nPix);
            dqPix = obj.experiment.dqPix;
            % Basis vectors for decomposing dQ: [ 3 x nPix ]
            qBasis_rad = normalizeVectorList(obj.qIMG_list);
            qBasis_par = normalizeVectorList(k0_mesh+0.5*obj.qIMG_list);
            qBasis_rot = cross(qBasis_par,qBasis_rad,2);

            % sigma values for describing peak shapes, dim: [nPix x 1]
            sigmaRad2rd = (obj.model.thetaRadial*obj.modqIMG_list).^2 + dqPix^2;
            sigmaPar2rd = (obj.model.thetaParallel*obj.modqIMG_list).^2 + ...
                          (obj.experiment.beamDivSigma*obj.modqIMG_list).^2 + obj.model.sigmaEwald^2;
            sigmaRot2rd = (obj.model.thetaRotation*obj.modqIMG_list).^2 + dqPix^2;
            

            invSigmaRad2rd = 0.5./sigmaRad2rd;
            invSigmaPar2rd = 0.5./sigmaPar2rd;
            invSigmaRot2rd = 0.5./sigmaRot2rd;
            
            disp(['nHKL = ',num2str(obj.nHKL)]);
  
            Vpixels = 1:obj.nPix;
            %sigmaPar is not included in normalization since detection integrates this
            %dimension
            Normalization = 1./sqrt(sigmaRad2rd.*sigmaRot2rd.*sigmaPar2rd);
            

            tic
            for indx = 1:obj.nHKL
                qHKL_mesh = meshgrid(obj.qHKL_list(indx,:),Vpixels);
                dQ = qHKL_mesh - obj.qIMG_list;
                dQrad2rd = dot(dQ,qBasis_rad,2).^2;
                dQpar2rd = dot(dQ,qBasis_par,2).^2;
                dQrot2rd = dot(dQ,qBasis_rot,2).^2;
                argExp = dQrad2rd.*invSigmaRad2rd + dQpar2rd.*invSigmaPar2rd + dQrot2rd.*invSigmaRot2rd;
                obj.WTS(:,indx) = Normalization.*exp(-argExp);
            end
            toc
        end %WTS
        

        
        function Fhkl = calcFhkls(obj,UVW,Uiso)
            [Umesh,modQmesh] = meshgrid(Uiso,obj.modQ_hkl);
            Fhkl_uvw = obj.FatomicT.*exp(-2*pi*(Umesh.*modQmesh).^2).*cos(2*pi*obj.HKL_list*UVW');
            Fhkl = sum(Fhkl_uvw,2);
        end
        
        function obj = calcIhkl_GS(obj)
            obj.Fhkl_GS = obj.calcFhkls(obj.model.UVW_GS, obj.model.Uiso_GS);
            obj.Ihkl_GS = abs(obj.Fhkl_GS).^2;
        end
        
        function obj = calcIhkl_ES(obj)
            obj.Fhkl_ES = obj.calcFhkls(obj.UVW_ES,obj.Uiso_ES);
            obj.Ihkl_ES = abs(obj.xExcited*obj.Fhkl_ES + (1-obj.xExcited)*obj.Fhkl_GS).^2;
            obj.Ihkl_diff = obj.Ihkl_ES - obj.Ihkl_GS;
        end
        
        function IMG = calcIMG(obj, Ihkls)
            IMG_1D = obj.WTS*Ihkls;
            IMG = imgaussfilt(reshape(IMG_1D, [obj.experiment.imgsize.y, obj.experiment.imgsize.x]), obj.experiment.spotSizePixSigma);
            if obj.model.reflectImage
                IMG = fliplr(IMG);
            end  
        end
        
        function obj = calcIMG_GS(obj)
            obj.IMG_GS = obj.calcIMG(obj.Ihkl_GS);
        end
        
        function obj = fitResults2diffPattern(obj,RESULTS,Tindex)
            obj.UVW_ES = RESULTS.UVW_T(:,:,Tindex);
            obj.Uiso_ES = RESULTS.Uiso_T(:,Tindex);
            obj.xExcited = RESULTS.xExcited;
            obj = calcIMG_diff(obj);
        end
        
        function obj = calcIMG_diff(obj)
            obj = obj.calcIhkl_ES();
            obj.IMG_diff = obj.calcIMG(obj.Ihkl_diff);
        end
        
        function obj = runSim_GS(obj)
            if isempty(obj.WTS)
                obj = obj.calculateWTS();
            end
            if isempty(obj.Ihkl_GS)
                obj = obj.calcIhkl_GS();
            end
            obj = calcIMG_GS(obj);
        end
        
        function I000 = findI000(obj)
            rows = 1:obj.nHKL;
            eq0 = obj.HKL_list == 0;
            indxHKL000 = rows(sum(eq0,2) == 3);
            I000 = obj.Ihkl_GS(indxHKL000);
        end
        
    end %methods
    
end %classdef
    
function HKL_list = createHKL_list(minhkls,maxhkls)
    Hvals = minhkls(1):maxhkls(1);
    Kvals = minhkls(2):maxhkls(2);
    Lvals = minhkls(3):maxhkls(3);
    [KKK,HHH,LLL] = meshgrid(Kvals,Hvals,Lvals); %note meshgrid takes order xyz=213
    nhkl = maxhkls-minhkls+1;
    nVals = nhkl(1)*nhkl(2)*nhkl(3);
    Hvec = reshape(HHH,[nVals,1]);
    Kvec = reshape(KKK,[nVals,1]);
    Lvec = reshape(LLL,[nVals,1]);
    HKL_list = [Hvec,Kvec,Lvec];
end


function Tmat = transformMatFrom_k0(k0)
    V3 = k0/norm(k0); % e-beam in z-direction
    Vy = cross(V3,[1,0,0]);
    V2 = Vy/norm(Vy);
    V1 = cross(V2,V3);
    Tmat = [V1;V2;V3];
end

function VLnorm = normalizeVectorList(VL)
    modVL = sqrt(sum(VL.^2,2));
    modVLmesh = meshgrid(modVL,1:3)';
    VLnorm = VL./modVLmesh;
end 