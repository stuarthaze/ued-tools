classdef DiffractionCalculator6
% Major revision compared to v5
% Can use arbitrary texturing axis
% Cartesians defined as:
% x || a
% y in ab plane

    
    properties
        crystal
        experiment
        % Unit cell properties     
        a
        b
        c
        alpha
        beta
        gamma
        alpha_deg
        beta_deg
        gamma_deg
        V
        axes
        T_FracCart     % xyz = uvw*T
        T_CartFrac
        T_hkl2xyz
        T_xyz2hkl
        T_lab2xyz
        T_xyz2lab
        T_tex2xyz
        T_xyz2tex
        % Structure
        nAt
        atomTypes
        atomicNumbers
        UVW
        XYZ
        occup
        Uiso
        Uijn
        % --- Reciprocal Space
        % All HKLs within resolution cut-off
        qMax
        nHKL_all
        HKL_all
        qHKL_all    % Cartesian
        modqHKL_all
        kHKL_all    % Cartesian
        khkl_all    % polar from beam direction
        usedHKL     % Logical
        % Reduced list
        nHKL
        HKL
        qHKL        % Cartesian Crystal frame
        qhkl        % Polar from beam/texture axis
        khkl        % Polar from beam direction
        modqHKL
        modqHKL2rd
        Fatomic     % 2d array [nAt x nHKL]
        FatomicT
        Fhkl
        Ihkl
        indxHKL000
        xExcited
        fAtomicType
        fAtomicComplex
        thetaParallel
        thetaRadial
        thetaRotation
        imgRotation
        dqPix
        dqMin
        dqMinParallel
        beamDivergence    % 1*sigma
%         UVW_GS
%         UVW_ES
%         Uiso_GS
%         Uiso_ES
%         Fhkl_GS
%         Fhkl_ES
%         Ihkl_GS
%         Ihkl_ES
%         Ihkl_diff
        % Experiment parameters
        nPix
        nx
        ny
        cenx
        ceny
        reflectImage
        k0_uvw
        k0_xyz  %Normalized
        modk0
        wavelength
        useTextureAxis
        textureAxis_xyz   % Axis for texture in Cartesian
        textureAxis_uvw
   
        qIMG       % Cartesian, Lab coordinate frame: Beam in [0,0,1] direction
        modqIMG
        qIMG_list
        qpix           % Polar wrt beam/texture axis
%         kpix           % Polar wrt beam
        maxqIMG
        modqIMG_list
        
        WTS %Matrix to weight reflections: [nPix x nHKL]
        IMG
        IMG_diff
        
    end %properties
    %---------------------------------------------------------------------
    
    methods
        % Constructor
        function obj = DiffractionCalculator6()
        end
        
        function obj = loadCrystalDataFromStruct(obj,CRYSTAL)
            obj.crystal = CRYSTAL;
            obj.a = CRYSTAL.a;
            obj.b = CRYSTAL.b;
            obj.c = CRYSTAL.c;
            
            obj.alpha = CRYSTAL.alpha;
            obj.beta = CRYSTAL.beta;
            obj.gamma = CRYSTAL.gamma;
            obj.alpha_deg = CRYSTAL.alpha_deg;
            obj.beta_deg = CRYSTAL.beta_deg;
            obj.gamma_deg = CRYSTAL.gamma_deg;
            
            obj.V = CRYSTAL.V;
            obj.axes = CRYSTAL.axes;
            obj.T_FracCart = CRYSTAL.T_FracCart;
            obj.T_CartFrac = CRYSTAL.T_CartFrac;
            obj.T_hkl2xyz = CRYSTAL.T_hkl2xyz;
            obj.T_xyz2hkl = CRYSTAL.T_xyz2hkl;
            
            obj.nAt = CRYSTAL.nAt;
            obj.UVW = CRYSTAL.UVW;
            obj.Uiso = CRYSTAL.Uiso;
            obj.atomTypes = CRYSTAL.atomTypes;
            obj.atomicNumbers = CRYSTAL.atomicNumbers;
        end
        
        function obj = readExperimentParameters(obj,EXPERIMENT)
            obj.experiment = EXPERIMENT;
            obj.wavelength = EXPERIMENT.wavelength;
            obj.modk0 = 1/obj.wavelength;
            obj.k0_uvw = EXPERIMENT.beamDir_uvw;
            obj.cenx = EXPERIMENT.cen.x;
            obj.ceny = EXPERIMENT.cen.y;
            obj.dqPix = EXPERIMENT.pixelSize*obj.modk0/(EXPERIMENT.magnification*EXPERIMENT.detectorDist);
            obj.dqMin = obj.dqPix;
            obj.beamDivergence = EXPERIMENT.beamDivergenceFWHM/2.35;
            obj = obj.calck0fromFractionals();
            obj.nx = EXPERIMENT.imgsize.x;
            obj.ny = EXPERIMENT.imgsize.y;
            obj.reflectImage = EXPERIMENT.reflectImage;
            obj.imgRotation = EXPERIMENT.imgRotation;
        end

        function obj = calck0fromFractionals(obj)
            k0 = obj.k0_uvw*obj.T_FracCart;
            obj.k0_xyz = k0/norm(k0);
        end
        
        function obj = calc_qIMG_labframe(obj)
            obj.nx = obj.nx;
            obj.ny = obj.ny;
            obj.nPix = obj.nx*obj.ny;
            scale = obj.dqPix;
            [XX,YY] = meshgrid(1:obj.nx,1:obj.ny);
            ZZ = obj.modk0*ones(obj.ny,obj.nx);
            dXX = scale*(XX-obj.cenx);
            dYY = scale*(YY-obj.ceny);
            if obj.reflectImage
                dXX = -dXX;
            end
            [el,az,r] = cart2sph(dXX,dYY,ZZ);
            r(:,:) = obj.modk0;
            [kx,ky,kz] = sph2cart(el,az,r);
            qXYZ = zeros([obj.ny,obj.nx,3]);
            qXYZ(:,:,1) = kx;
            qXYZ(:,:,2) = ky;
            qXYZ(:,:,3) = kz-obj.modk0;
            obj.qIMG = qXYZ;
            obj = obj.calc_modqIMG();
        end
        
        function obj = calc_modqIMG(obj)
            obj.modqIMG = sqrt(sum(obj.qIMG.^2,3));
            obj.qMax = max(max(obj.modqIMG));
        end
        
        
        function obj = generateHKLlist(obj)
            HKLlist = createHKL_list_v2(obj.T_xyz2hkl,obj.qMax);
            qHKLlist = HKLlist*obj.T_hkl2xyz;
            modqlist = sqrt(sum(qHKLlist.^2,2));
            indx2use = modqlist <= obj.qMax;
            obj.HKL_all = HKLlist(indx2use,:);
            obj.modqHKL_all = modqlist(indx2use);
            obj.qHKL_all = obj.HKL_all*obj.T_hkl2xyz;
            obj.kHKL_all = obj.qHKL_all*obj.T_xyz2lab;
            obj.kHKL_all(:,3) = obj.kHKL_all(:,3) + obj.modk0;
            obj.khkl_all = cart2sphPolar(obj.kHKL_all);
            obj.nHKL_all = size(obj.HKL_all,1);
        end
        
        function modkHKL = hkl2modk(obj,hkl)
            qvals = hkl*obj.T_hkl2xyz*obj.T_xyz2lab;
            kvals = qvals;
            kvals(:,3) = qvals(:,3)+ obj.modk0;
            k_polar = cart2sphPolar(kvals);
            modkHKL = k_polar.r;
        end
                
        function obj = cropHKLs(obj)
            if (obj.thetaParallel+obj.thetaRotation) > 0.5*pi
                obj.usedHKL = (1:obj.nHKL_all) > 0;
            else
                cutoff = 2*obj.qMax*sqrt(obj.thetaParallel^2 + obj.thetaRotation^2)+ 3*obj.dqMin;
                dk_hkl = obj.khkl_all.r - obj.modk0;
                obj.usedHKL = (dk_hkl >= -cutoff) & (dk_hkl <= cutoff);
            end
            obj = obj.selectHKLs();
        end
        
        function obj = selectHKLs(obj)
            indices = obj.usedHKL;
            obj.nHKL = sum(indices);
            obj.HKL = obj.HKL_all(indices,:);
            obj.qHKL = obj.qHKL_all(indices,:);
            obj.modqHKL = obj.modqHKL_all(indices,:);
            obj.modqHKL2rd = obj.modqHKL.^2;
        end
        
        % ---- Transformation matrices ------------
        function obj = calc_T_lab2xyz(obj,vx)
            obj.T_lab2xyz = calculateTransformationMatrix(obj.k0_xyz,vx);
            obj.T_xyz2lab = eye(3)/obj.T_lab2xyz;
        end
        
        function obj = calculateTranformationForTextureAxis(obj,vx)
            obj.textureAxis_xyz = obj.textureAxis_uvw*obj.T_FracCart;
            if nargin < 2
                obj.T_tex2xyz = calculateTransformationMatrix(obj.textureAxis_xyz);
            else
                obj.T_tex2xyz = calculateTransformationMatrix(obj.textureAxis_xyz,vx);
            end
            obj.T_xyz2tex = eye(3)/obj.T_tex2xyz;
        end
        % -------
        function obj = calc_qPolar(obj)
            obj.qIMG_list = reshape(obj.qIMG,[obj.nPix,3]);
            if obj.useTextureAxis
                obj = obj.calculateTranformationForTextureAxis();
                obj.qpix = cart2sphPolar(obj.qIMG_list*obj.T_lab2xyz*obj.T_xyz2tex);
                obj.qhkl = cart2sphPolar(obj.qHKL*obj.T_xyz2tex);
            else
                %Using Lab-frame
                obj.qpix = cart2sphPolar(obj.qIMG_list);
                obj.qhkl = cart2sphPolar(obj.qHKL*obj.T_xyz2lab);
            end
        end
        
        function obj = calculateFatomic(obj)
            if strcmp(obj.fAtomicType,'Kirkland')
                ScatFactTable = load('KirklandTable','-ascii');
                modS = 2*pi*obj.modqHKL;
                obj.Fatomic = zeros([obj.nAt,obj.nHKL]);
                for s = 1:obj.nHKL
                    for atm = 1:obj.nAt
                        obj.Fatomic(atm,s) = KirklandTable2f(obj.atomicNumbers(atm),modS(s),ScatFactTable);
                    end
                end
            elseif strcmp(obj.fAtomicType,'EDICO')
                modS_hkl = modQ_hkl*2*pi;
                [F_mod,F_phi] = interp_f_from_tables_EDICO(CRYSTAL.atomicNumbers,modS_hkl,EXPERIMENT.voltage_kV);
                if MODEL.fAtomicComplex
                    obj.Fatomic = (F_mod.*exp(1i*F_phi))';
                else
                    obj.Fatomic = F_mod';
                end
            elseif isempty(obj.fAtomicType)
                disp('Need to specify atomic scattering data source');
            end
            obj.FatomicT = obj.Fatomic';
        end
        
        
        

        


















%             obj.model = MODEL;
%  
%             % Calculate detector geometry in terms of q
%             imgsize = EXPERIMENT.imgsize;
%             cen = EXPERIMENT.cen;
%             pixelscale = EXPERIMENT.pixelSize/(EXPERIMENT.detectorDist*EXPERIMENT.magnification); % angle between pixels at origin
%             
%             modk0 = 1/EXPERIMENT.wavelength;
%             obj.experiment.dqPix = modk0*pixelscale;
%             beamSizeSigma = EXPERIMENT.beamSizeFWHM/2.35;
%             beamDivSigma = beamSizeSigma/(EXPERIMENT.detectorDist*EXPERIMENT.magnification);
%             obj.experiment.beamDivSigma = beamDivSigma;
%             obj.experiment.spotSizePixSigma = beamSizeSigma/EXPERIMENT.pixelSize;
% 
% %             [qIMG, modQ] = calcQimg(imgsize.x,imgsize.y,pixelscale,cen,modk0);
%             qIMG = calcqIMG_v2(imgsize.x,imgsize.y,pixelscale,cen,modk0,EXPERIMENT.reflectImage);
%             
%             nPix = imgsize.x*imgsize.y;
%             qIMG_L = reshape(qIMG,[nPix,3]); % q-vectors for img
%             qpix = cart2sphPolar(qIMG_L);
% 
% %             modqIMG_L = reshape(modQ,[nPix,1]);
%             maxqIMG = max(qpix.r);
%             kpix = qIMG_L + [0,0,modk0];
%             
%             %Calculate transformation matrices
%             [T_h2q,T_q2h] = calcTransformMatrix_h2q(obj.experiment.beamDir_uvw,obj.crystal);
% 
%             HKL_maxList = createHKL_list_v2(CRYSTAL.T_xyz2hkl,maxqIMG);
%             nHKL_maxList = size(HKL_maxList,1)
%             qHKL_maxList = HKL_maxList*T_h2q;  
%             kHKL_maxList = qHKL_maxList + meshgrid([0,0,modk0],1:nHKL_maxList);
%             mod_kHKL_maxList = sqrt(sum(kHKL_maxList.^2,2));
%             dmod_kHKL_maxList = mod_kHKL_maxList - modk0;
%             
%             nSigmaCutoff = 2;
%             maxSigmaParallel = sqrt((MODEL.thetaParallel*maxqIMG)^2 + (beamDivSigma*maxqIMG)^2 + MODEL.sigmaEwald^2);
%             dkCutoff = nSigmaCutoff*maxSigmaParallel;
%             modqHKL_maxList = sqrt(sum(qHKL_maxList.^2,2));
% 
%             indices2keep = (dmod_kHKL_maxList >= -dkCutoff) & (dmod_kHKL_maxList <= dkCutoff) & (modqHKL_maxList <= maxqIMG);
%             HKLs = HKL_maxList(indices2keep,:);
%             qHKLs = qHKL_maxList(indices2keep,:);
%             khkl.list = kHKL_maxList(indices2keep,:);
%             khkl.dk = dmod_kHKL_maxList(indices2keep,:);
%             modQ2rd_hkl = sum(qHKLs.^2,2);
%             modQ_hkl = sqrt(modQ2rd_hkl);
%             
%             %Convert q to spherical polar coords
%             qhkl = cart2sphPolar(qHKLs);
%             
%             % Calculate Fatomic for HKLs
%             if strcmp(MODEL.fAtomicType,'Kirkland')
%                 Fatomic = calculateAtomicF_Kirkland_v2(CRYSTAL,HKLs);
%             elseif strcmp(MODEL.fAtomicType,'EDICO')
%                 modS_hkl = modQ_hkl*2*pi;
%                 [F_mod,F_phi] = interp_f_from_tables_EDICO(CRYSTAL.atomicNumbers,modS_hkl,EXPERIMENT.voltage_kV);
%                 if MODEL.fAtomicComplex
%                     Fatomic = (F_mod.*exp(1i*F_phi))';
%                 else
%                     Fatomic = F_mod';
%                 end
%             elseif ~exist('MODEL.fAtomicType','var')
%                 disp('Need to specify atomic scattering data source');
%             end 
% 
%             %Set object properties required later
%             obj.HKLs = HKLs;
%             obj.nHKL = size(HKLs,1);
%             obj.qHKLs = qHKLs;
%             obj.qhkl = qhkl;
%             obj.khkl = khkl;
%             obj.qpix = qpix;
%             obj.kpix = kpix;
%             obj.modQ_hkl = modQ_hkl;
%             obj.modQ2rd_hkl = modQ2rd_hkl;
%             obj.k0 = [0,0,modk0];
%             obj.modk0 = modk0;
%             obj.nPix = nPix;
%             obj.imgsize = obj.experiment.imgsize;
%             obj.qIMG = qIMG;
%             obj.qIMG_list = qIMG_L;
%             obj.maxqIMG = maxqIMG;
% %             obj.modqIMG_list = modqIMG_L;
%             obj.Fatomic = Fatomic;
%             obj.FatomicT = Fatomic';
% 
%         end %constructor
%         
        %------------------------------------------------------------------
        
        function obj = calculateWTS_polar(obj,use000)
            obj = obj.findIndxHKL000();
            obj.WTS = zeros([obj.nPix,obj.nHKL]);
            dqMin2rd = obj.dqMin^2;
            if isempty(obj.dqMinParallel)
                dqMinPar = obj.dqMin;
            else
                dqMinPar = obj.dqMinParallel;
            end
            dqMinPar2rd = dqMinPar^2;
            
            
            for indx = 1:obj.nHKL
                modq = obj.qhkl.r(indx);
                cosElev = cos(obj.qhkl.el(indx));
                sig_qRad_2rd = (obj.thetaRadial*modq)^2 + 0.5*dqMin2rd;
                sig_qRot_2rd = (obj.thetaRotation*modq*cosElev)^2 + 0.5*dqMin2rd;
                sig_qPar_2rd = (obj.thetaParallel*modq)^2 + (obj.beamDivergence*modq)^2 + 0.5*dqMinPar2rd;
%                     Normalization = 1/sqrt(sig_qRad_2rd*sig_qRot_2rd*sig_qPar_2rd);

                % calc dq for all pixels and HKL(indx)
                dq_rad = obj.qpix.r - obj.qhkl.r(indx);
                dq_az  = obj.qpix.az- obj.qhkl.az(indx)+obj.imgRotation; %Doesn't work with Texture axis + obj.model.imgRotation; %Note: Could also add an image rotation here
                dq_par = (obj.qpix.el- obj.qhkl.el(indx))*modq;
%                     %Correct az values for errors of 2*pi
%                 dq_rot = mod(abs(dq_az),2*pi)*modq;
                dq_az = mod(abs(dq_az),2*pi);
                dq_az_neg = 2*pi - dq_az;
                dq_rot = min(dq_az,dq_az_neg)*modq*cosElev;

                if ~(indx == obj.indxHKL000)
                    argExp = (dq_rad.^2*0.5/sig_qRad_2rd)+(dq_rot.^2*0.5/sig_qRot_2rd)+(dq_par.^2*0.5/sig_qPar_2rd);
                    Normalization = dqMin2rd/sqrt(sig_qRad_2rd*sig_qRot_2rd);
                elseif use000
                    argExp = dq_rad.^2*0.5/sig_qRad_2rd;
                    Normalization = dqMin2rd/sig_qRad_2rd;
                else
                    Normalization = 0;
                end
                obj.WTS(:,indx) = Normalization*exp(-argExp);
            end
        end %WTS polar
        %------------------------------------------------------------------
        
        function Fhkl = calcFhkls(obj,UVW,Uiso)
            if isempty(Uiso)
                Fhkl_uvw = obj.FatomicT.*cos(2*pi*obj.HKL*UVW');
            else
                [Umesh,modQ2rdMesh] = meshgrid(Uiso,obj.modqHKL2rd);
                Fhkl_uvw = obj.FatomicT.*exp(-2*pi^2*Umesh.*modQ2rdMesh).*cos(2*pi*obj.HKL*UVW');
            end
            Fhkl = sum(Fhkl_uvw,2);
        end
        
        function obj = calcIhkl_GS(obj)
            obj.Fhkl = obj.calcFhkls(obj.UVW, obj.Uiso);
            obj.Ihkl = abs(obj.Fhkl).^2;
        end
        
        function obj = calcIhkl_ES(obj)
            obj.Fhkl_ES = obj.calcFhkls(obj.UVW_ES,obj.Uiso_ES);
            obj.Ihkl_ES = abs(obj.xExcited*obj.Fhkl_ES + (1-obj.xExcited)*obj.Fhkl_GS).^2;
            obj.Ihkl_diff = obj.Ihkl_ES - obj.Ihkl_GS;
        end
        
        function IMG = calcIMG(obj, Ihkls)
            IMG_1D = obj.WTS*Ihkls;
            IMG = reshape(IMG_1D, [obj.experiment.imgsize.y, obj.experiment.imgsize.x]);
        end
        
        function obj = calcIMG_GS(obj)
            obj.IMG = obj.calcIMG(obj.Ihkl);
        end
              
        function obj = calcIMG_diff(obj)
            obj = obj.calcIhkl_ES();
            obj.IMG_diff = obj.calcIMG(obj.Ihkl_diff);
        end
        
        function obj = fitResults2diffPattern(obj,RESULTS,Tindex)
            obj.UVW_ES = RESULTS.UVW_T(:,:,Tindex);
            obj.Uiso_ES = RESULTS.Uiso_T(:,Tindex);
            obj.xExcited = RESULTS.xExcited;
            obj = calcIMG_diff(obj);
        end

        %-------------------------------------------------------------
        % Run Simulation
        function obj = runSim_GS(obj,options)
            if options.timeSimulation
                tic
            end
            
            if isempty(obj.WTS) || options.recalculateWTS
                obj = obj.calculateWTS_polar(options.use000);
            end
            
            if isempty(obj.UVW)
                disp('Enter fractional coordinates');
            end
            
            if isempty(obj.Uiso)
                obj.Uiso = zeros([obj.nAt,1]);
            end
            
            obj = obj.calcIhkl_GS();

            obj = obj.calcIMG_GS();
            
            if options.beamSizeSigmaPix
                obj.IMG = imgaussfilt(obj.IMG,options.beamSizeSigmaPix);
            end

            if options.timeSimulation
                toc
            end
            if options.displayIMG
                obj.displayIMG(options.scaleIMG)
            end
        end
        %----------------------------------------------------------
        function I000 = returnI000(obj)
            obj = obj.findIndxHKL000();
            I000 = obj.Ihkl(obj.indxHKL000);
        end
        
        function obj = findIndxHKL000(obj)
            rows = 1:obj.nHKL;
            eq0 = obj.HKL == 0;
            obj.indxHKL000 = rows(sum(eq0,2) == 3);
        end
        
        function result = findHKLrow(obj,hkl)
            rows = 1:obj.nHKL;
            hklmatch = obj.HKL == hkl;
            result = rows(sum(hklmatch,2) == 3);
        end
        
        function result = getIhklGS(obj,hkl)
            indxHKL = obj.findHKLrow(hkl);
            result = obj.Ihkl_GS(indxHKL);
        end
        function result = get_Ihkl_diff(obj,hkl)
            indxHKL = obj.findHKLrow(hkl);
            result = obj.Ihkl_diff(indxHKL);
        end
        function result = getRatioIhkl(obj,hkl)
            indxHKL = obj.findHKLrow(hkl);
            Igs = obj.Ihkl_GS(indxHKL);
            Ies = obj.Ihkl_ES(indxHKL);
            result = Ies/Igs;
        end
        
        function obj = readSimModel(obj,MODEL)
            obj.model = MODEL;
        end
        
        function displayIMG(obj,scale)
            if (nargin < 2) || ~scale
                scale = max(max(obj.IMG));
            end
            figure();
            imagesc(obj.IMG,[0,scale]); 
            pbaspect([obj.experiment.imgsize.x,obj.experiment.imgsize.y,1]); 
            colormap(flipud(gray));
        end
        
        function displaySqrtIMG(obj,scale)
            figure();
            sqrtIMG = sqrt(obj.IMG);
            if nargin < 2
                scale = max(max(sqrtIMG));
            end
            imagesc(sqrtIMG,[0,scale]); 
            pbaspect([obj.experiment.imgsize.x,obj.experiment.imgsize.y,1]); 
            colormap(flipud(gray));
        end
        
        function displayDifferenceIMG(obj,scale)
            if isempty(obj.IMG_diff)
                disp('No IMG_diff');
                return
            elseif (nargin < 2) || ~scale
                scale = 0.5*(max(max(obj.IMG_diff))-min(min(obj.IMG_diff)));
            end
            figure();
            imagesc(obj.IMG_diff,[-scale,scale]); 
            pbaspect([obj.imgsize.x,obj.imgsize.y,1]); 
            colormap('gray');            
        end
        
        function pix = getPixelDetails(obj,X,Y,nCont)
            if nargin == 3
                nCont = 1;
            end
            pix.x = X;
            pix.y = Y;
            pix.Igs = obj.IMG(Y,X);
            pix.indx = obj.xy2indx(X,Y);
            Wpix = obj.WTS(pix.indx,:);
            [Wsorted, pix.indxHKL_sorted] = sort(Wpix,'descend');
            pix.indxHKL_max = pix.indxHKL_sorted(1);
            pix.WtsMax = Wsorted(1:nCont);
            pix.HKL = obj.HKL(pix.indxHKL_sorted(1:nCont),:);
            pix.Fhkl_maxWts = obj.Fhkl(pix.indxHKL_sorted(1:nCont))';
            pix.q = obj.qIMG_list(pix.indx,:);
            pix.q_r = obj.qpix.r(pix.indx);
%             if obj.model.reflectImage
%                 disp('Warning: Image reflected');
%             end
        end
        
        function indx = xy2indx(obj,X,Y)
            indx = Y+((X-1)*obj.ny);
        end
        
        function [X,Y] = indx2xy(obj,indx)
            X = ceil(indx/obj.ny);
            Y = mod(indx-1,obj.ny) + 1;
        end

            
    end %methods
    
end %classdef

%-------------------------------------------------------
% Utility Functions
% %-------------------------------------------------------
% function [Q,Qmod] = calcQimg(nx,ny,pixScale,cen,modk0,reflectImage)
%     Q = zeros([ny,nx,3]);
%     [XX,YY] = meshgrid(1:nx,1:ny);
%     dXX = XX-cen.x;
%     if reflectImage
%         dXX = -dXX;
%     end  
%     dYY = YY-cen.y;
%     RR = hypot(dXX,dYY);
%     AA = atan(pixScale*RR); %Scattering angles
%     Qmod = 2*modk0*sin(0.5*AA);
%     Q(:,:,3) = -Qmod.*sin(0.5*AA);
%     qRR = Qmod.*cos(0.5*AA);
%     qXX = qRR.*dXX./RR;
%     qYY = qRR.*dYY./RR;
%     Q(:,:,1) = qXX;
%     Q(:,:,2) = qYY;
% end

function qXYZ = calcqIMG_v2(nx,ny,scale,cen,modk0,reflectImage)
    [XX,YY] = meshgrid(1:nx,1:ny);
    ZZ = ones(ny,nx);
    dXX = scale*(XX-cen.x);
    dYY = scale*(YY-cen.y);
    if reflectImage
        dXX = -dXX;
    end
    [el,az,r] = cart2sph(dXX,dYY,ZZ);
    r(:,:) = modk0;
    [kx,ky,kz] = sph2cart(el,az,r);
    qXYZ = zeros([ny,nx,3]);
    qXYZ(:,:,1) = kx;
    qXYZ(:,:,2) = ky;
    qXYZ(:,:,3) = kz-modk0;
end

function HKL_list = createHKL_list_v2(T_xyz2hkl,resolution)
    maxhkls = round(resolution*[1,1,1]*T_xyz2hkl);
    Hvals = -maxhkls(1):maxhkls(1);
    Kvals = -maxhkls(2):maxhkls(2);
    Lvals = -maxhkls(3):maxhkls(3);
    [KKK,HHH,LLL] = meshgrid(Kvals,Hvals,Lvals); %note meshgrid takes order xyz=213
    nhkl = (2*maxhkls) + [1,1,1];
    nVals = nhkl(1)*nhkl(2)*nhkl(3);
    Hvec = reshape(HHH,[nVals,1]);
    Kvec = reshape(KKK,[nVals,1]);
    Lvec = reshape(LLL,[nVals,1]);
    HKL_list = [Hvec,Kvec,Lvec];
end
% 
% function VLnorm = normalizeVectorList(VL)
%     modVL = sqrt(sum(VL.^2,2));
%     modVLmesh = meshgrid(modVL,1:3)';
%     VLnorm = VL./modVLmesh;
% end 

function R = cart2sphPolar(XYZ)
    % Input: XYZ is an array [nAtom x 3]
    % OUTPUT:
    % R.r = radius
    % R.az = azimuthal is atan(y/x)
    % R.el = elevation is angle from xy plane = pi/2 - theta
    [R.az,R.el,R.r] = cart2sph(XYZ(:,1),XYZ(:,2),XYZ(:,3));
end

function XYZ = sphPolar2cart(R)
    [x,y,z] = sph2cart(R.az,R.el,R.r);
    XYZ = [x,y,z];
end

function Tmat = calculateTransformationMatrix(V3z,V1x)
% Input: V3z - new z-axis in original basis
%        V1x - Vector to define new x-axis perpendicular to V3
    V3 = V3z/norm(V3z); % e-beam in z-direction
    if nargin == 2
        Vxplane = V1x;
        Vy = cross(V3,Vxplane);
    else
        Vy = cross(V3,[1,0,0]);
        if norm(Vy) < 0.5
            Vy = cross(V3,[0,1,0]);
        end
    end    
    if norm(Vy) < 0.5
        Vy = cross(V3,[0,1,0]);
    end
    V2 = Vy/norm(Vy);
    V1 = cross(V2,V3);
    Tmat = [V1;V2;V3];
end

function Tmat = calculateTransformMatrix_Lab2CrystalCart(k0)
% Input: k0 beam direction in Cartesian coords (crystal frame)
    V3 = k0/norm(k0); % e-beam in z-direction
    Vy = cross(V3,[1,0,0]);
    if norm(Vy) < 0.5
        Vy = cross(V3,[0,1,0]);
    end
    V2 = Vy/norm(Vy);
    V1 = cross(V2,V3);
    Tmat = [V1;V2;V3];
end

function [Th2q,Tq2h] = calcTransformMatrix_h2q(beamDir_uvw,crystal)
    beamDir_xyz = beamDir_uvw*crystal.T_FracCart;
    T_lab2cryst = calculateTransformMatrix_Lab2CrystalCart(beamDir_xyz);
    Tq2h = T_lab2cryst*crystal.T_xyz2hkl;
    Th2q = eye(3)/Tq2h;
end