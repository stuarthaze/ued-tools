classdef DiffractionCalculator2
% Uses lab-frame coordinates where the electron beam,
% k0 is in the z-direction, (0,0,1) 

    
    properties
        crystal
        experiment
        model
        %Calculated internally
        nHKL
        HKLs
        qHKLs
        qhkl
        khkl
        modQ_hkl
        Fatomic
        FatomicT
        indxHKL000
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
        imgsize
        k0
        modk0
        qIMG
        qIMG_list
        qpix
        kpix
        maxqIMG
        modqIMG_list
        
        WTS %Matrix to weight reflections: [nPix x nHKL]
        IMG
        IMG_diff
        
    end %properties
    %---------------------------------------------------------------------
    
    methods
        % Constructor
        function obj = DiffractionCalculator2(CRYSTAL,EXPERIMENT,MODEL)
            obj.crystal = CRYSTAL;
            obj.experiment = EXPERIMENT;
            obj.model = MODEL;
            
            %--------------------------------------------------------
            %            Nested functions for constructor
            
            function [Q,Qmod] = calcQimg()
                phi = EXPERIMENT.imgRotation;
                Q = zeros([imgsize.y,imgsize.x,3]);
                [XX,YY] = meshgrid(1:imgsize.x,1:imgsize.y);
                dXX = XX-cen.x;
                dYY = YY-cen.y;
                RR = hypot(dXX,dYY);
                AA = atan(pixelscale*RR); %Scattering angles
                Qmod = 2*modk0*sin(0.5*AA);
                Q(:,:,3) = -Qmod.*sin(0.5*AA);
                qRR = Qmod.*cos(0.5*AA);
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

            %Calculate transformation matrix uvw -> q (lab frame)
            beamDir_xyz = EXPERIMENT.beamDir_uvw*CRYSTAL.T_FracCart;
            k0_xyz_norm = beamDir_xyz/norm(beamDir_xyz);
            k0_crystalFrame = k0_xyz_norm*modk0;
%             k0_normT = k0_xyz_norm';
            k0 = [0,0,modk0];
            [qIMG, modQ] = calcQimg();


            % Transformation Matrices (h,k,l) -> (qx,qy,qz)
            T_lab2cryst = calculateTransformMatrix_Lab2CrystalCart(k0_crystalFrame);
            T_q2h = T_lab2cryst*CRYSTAL.T_xyz2hkl;
            T_h2q = inv(T_q2h);
            
            % Calculate hkl values to use
            nPix = imgsize.x*imgsize.y;
            qIMG_L = reshape(qIMG,[nPix,3]); % q-vectors for img
            modqIMG_L = reshape(modQ,[nPix,1]);
            maxqIMG = max(modqIMG_L);
            kpix = qIMG_L + [0,0,modk0];
% 
%             HKL_IMG_L = round(qIMG_L*T_q2h);
%             HKL_minvals = min(HKL_IMG_L);
%             HKL_maxvals = max(HKL_IMG_L);
%             HKL_maxList = createHKL_list(HKL_minvals,HKL_maxvals);

            HKL_maxList = createHKL_list_v2(CRYSTAL.T_xyz2hkl,maxqIMG);
            nHKL_maxList = size(HKL_maxList,1)
            qHKL_maxList = HKL_maxList/T_q2h; 
            kHKL_maxList = qHKL_maxList + meshgrid([0,0,modk0],1:nHKL_maxList);
            mod_kHKL_maxList = sqrt(sum(kHKL_maxList.^2,2));
            dmod_kHKL_maxList = mod_kHKL_maxList - modk0;
            
%             qHKLs_dotk0 = qHKL_maxList*k0_normT;
%             qIMG_dotk0 = qIMG_L*k0_normT;
            nSigmaCutoff = 2;
            maxSigmaParallel = sqrt((MODEL.thetaParallel*maxqIMG)^2 + (beamDivSigma*maxqIMG)^2 + MODEL.sigmaEwald^2);
            dkCutoff = nSigmaCutoff*maxSigmaParallel;
%             qdotk0_min = min(qIMG_dotk0)-qCutoff
%             qdotk0_max = max(qIMG_dotk0)+qCutoff
             modqHKL_maxList = sqrt(sum(qHKL_maxList.^2,2));

            indices2keep = (dmod_kHKL_maxList >= -dkCutoff) & (dmod_kHKL_maxList <= dkCutoff) & (modqHKL_maxList <= maxqIMG);
            HKLs = HKL_maxList(indices2keep,:);
            qHKLs = qHKL_maxList(indices2keep,:);
            khkl.list = kHKL_maxList(indices2keep,:);
            khkl.dk = dmod_kHKL_maxList(indices2keep,:);
            modQ_hkl = sqrt(sum(qHKLs.^2,2));
            
            %Convert q to spherical polar coords
            qhkl = cart2sphPolar(qHKLs);
            qpix = cart2sphPolar(qIMG_L);
            
            % Calculate Fatomic for HKLs
            if strcmp(MODEL.fAtomicType,'Kirkland')
                Fatomic = calculateAtomicF_Kirkland(CRYSTAL,HKLs);
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
            obj.HKLs = HKLs;
            obj.nHKL = size(HKLs,1);
            obj.qHKLs = qHKLs;
            obj.qhkl = qhkl;
            obj.khkl = khkl;
            obj.qpix = qpix;
            obj.kpix = kpix;
            obj.modQ_hkl = modQ_hkl;
            obj.k0 = k0;
            obj.modk0 = modk0;
            obj.nPix = nPix;
            obj.imgsize = obj.experiment.imgsize;
            obj.qIMG = qIMG;
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
                qHKL_mesh = meshgrid(obj.qHKLs(indx,:),Vpixels);
                dQ = qHKL_mesh - obj.qIMG_list;
                dQrad2rd = dot(dQ,qBasis_rad,2).^2;
                dQpar2rd = dot(dQ,qBasis_par,2).^2;
                dQrot2rd = dot(dQ,qBasis_rot,2).^2;
                argExp = dQrad2rd.*invSigmaRad2rd + dQpar2rd.*invSigmaPar2rd + dQrot2rd.*invSigmaRot2rd;
                obj.WTS(:,indx) = Normalization.*exp(-argExp);
            end
            toc
        end %WTS
        %------------------------------------------------------------------
        
        function obj = calculateWTS_polar(obj)
            if isempty(obj.indxHKL000)
                obj = obj.findIndxHKL000();
            end
            if isempty(obj.WTS)
                obj.WTS = zeros([obj.nPix,obj.nHKL]);
            end
            dqPix = obj.experiment.dqPix;
            
            for indx = 1:obj.nHKL
                if ~(indx == obj.indxHKL000)
                    modq = obj.qhkl.r(indx);
                    sig_qRad_2rd = (obj.model.thetaRadial*modq)^2 + dqPix^2;
                    sig_qRot_2rd = (obj.model.thetaRotation*modq)^2 + dqPix^2;
                    sig_qPar_2rd = (obj.model.thetaParallel*modq)^2 + (obj.experiment.beamDivSigma*modq)^2 + obj.model.sigmaEwald^2;
                    Normalization = 1/sqrt(sig_qRad_2rd*sig_qRot_2rd*sig_qPar_2rd);
                    % calc dq for all pixels and HKL(indx)
                    dq_rad = obj.qpix.r - obj.qhkl.r(indx);
                    dq_az  = obj.qpix.az- obj.qhkl.az(indx) + obj.model.imgRotation; %Note: Could also add an image rotation here
                    dq_par = (obj.qpix.el- obj.qhkl.el(indx))*modq;
%                     %Correct az values for errors of 2*pi
                    dq_rot = mod(abs(dq_az),2*pi)*modq;
  
                    argExp = (dq_rad.^2*0.5/sig_qRad_2rd)+(dq_rot.^2*0.5/sig_qRot_2rd)+(dq_par.^2*0.5/sig_qPar_2rd);
                    obj.WTS(:,indx) = Normalization*exp(-argExp);
                end
            end
        end %WTS polar
        %------------------------------------------------------------------
        
        function Fhkl = calcFhkls(obj,UVW,Uiso)
            [Umesh,modQmesh] = meshgrid(Uiso,obj.modQ_hkl);
            Fhkl_uvw = obj.FatomicT.*exp(-2*pi*(Umesh.*modQmesh).^2).*cos(2*pi*obj.HKLs*UVW');
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
            obj.IMG = obj.calcIMG(obj.Ihkl_GS);
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

        function obj = runSim_GS(obj,options)
            if options.timeSimulation
                tic
            end
            
            if isempty(obj.WTS) || options.recalculateWTS
                obj = obj.calculateWTS_polar();
            end
            if isempty(obj.Ihkl_GS)
                obj = obj.calcIhkl_GS();
            end

            obj = calcIMG_GS(obj);

            if options.timeSimulation
                toc
            end
            if options.displayIMG
                obj.displayIMG(options.scaleIMG)
            end
        end
        
        function I000 = returnI000(obj)
            if isempty(obj.indxHKL000)
                obj.findIndxHKL000();
            end
            I000 = obj.Ihkl_GS(obj.indxHKL000);
        end
        
        function obj = findIndxHKL000(obj)
            rows = 1:obj.nHKL;
            eq0 = obj.HKLs == 0;
            obj.indxHKL000 = rows(sum(eq0,2) == 3);
        end
        
        function obj = readSimModel(obj,MODEL)
            obj.model = MODEL;
        end
        
        function displayIMG(obj,scale)
            if (nargin < 2) || ~scale
                scale = max(max(obj.IMG));
            end
            figure();
            imagesc(-obj.IMG,[-scale,0]); 
            pbaspect([obj.experiment.imgsize.x,obj.experiment.imgsize.y,1]); 
            colormap('gray');
        end
        
        function displaySqrtIMG(obj)
            figure();
            sqrtIMG = sqrt(obj.IMG);
            imagesc(-sqrtIMG,[-max(max(sqrtIMG)),0]); 
            pbaspect([obj.experiment.imgsize.x,obj.experiment.imgsize.y,1]); 
            colormap('gray');
        end
        
        function displayDifferenceIMG(obj,scale)
            if isempty(obj.IMG_diff)
                disp('No IMG_diff');
                return
            elseif (nargin < 2) || ~scale
                scale = 0.5*(max(max(obj.IMG_diff))-min(min(obj.IMG_diff)));
            end
            figure();
            imagesc(-obj.IMG_diff,[-scale,scale]); 
            pbaspect([obj.imgsize.x,obj.imgsize.y,1]); 
            colormap('gray');            
        end
        
        function pix = getPixelDetails(obj,X,Y,nCont)
            pix.x = X;
            pix.y = Y;
            pix.indx = obj.xy2indx(X,Y);
            Wpix = obj.WTS(pix.indx,:);
            [Wsorted, pix.indxHKL_sorted] = sort(Wpix,'descend');
            pix.indxHKL_max = pix.indxHKL_sorted(1);
            pix.Wmax = Wsorted(1:nCont);
            pix.HKLs = obj.HKLs(pix.indxHKL_sorted(1:nCont),:);
            pix.Fhkl_maxWts = obj.Fhkl_GS(pix.indxHKL_sorted(1:nCont))';
            pix.q = obj.qIMG_list(pix.indx,:);
            pix.q_r = obj.qpix.r(pix.indx);
            if obj.model.reflectImage
                disp('Warning: Image reflected');
            end
        end
        
        function indx = xy2indx(obj,X,Y)
            indx = Y+((X-1)*obj.imgsize.y);
        end
        
        function [X,Y] = indx2xy(obj,indx)
            [XX,YY] = meshgrid(1:obj.imgsize.x,1:obj.imgsize.y);
            X = XX(indx);
            Y = YY(indx);
        end
            
    end %methods
    
end %classdef

%-------------------------------------------------------
% Utility Functions
%-------------------------------------------------------
function HKL_list = createHKL_list(minhkls,maxhkls)
    Hvals = minhkls(1):maxhkls(1);
    Kvals = minhkls(2):maxhkls(2);
    Lvals = minhkls(3):maxhkls(3);
    [KKK,HHH,LLL] = meshgrid(Kvals,Hvals,Lvals); %note meshgrid takes order xyz=213
    nhkl = (maxhkls-minhkls) + [1,1,1];
    nVals = nhkl(1)*nhkl(2)*nhkl(3);
    Hvec = reshape(HHH,[nVals,1]);
    Kvec = reshape(KKK,[nVals,1]);
    Lvec = reshape(LLL,[nVals,1]);
    HKL_list = [Hvec,Kvec,Lvec];
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

function VLnorm = normalizeVectorList(VL)
    modVL = sqrt(sum(VL.^2,2));
    modVLmesh = meshgrid(modVL,1:3)';
    VLnorm = VL./modVLmesh;
end 

function Tmat = calculateTransformMatrix_Lab2CrystalCart(k0)
    V3 = k0/norm(k0); % e-beam in z-direction
    Vy = cross(V3,[1,0,0]);
    V2 = Vy/norm(Vy);
    V1 = cross(V2,V3);
    Tmat = [V1;V2;V3];
end

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
    XYZ = [x;y;z];
end
