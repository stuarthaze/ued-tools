function IMG = calcCrystalDiffractionPattern_v2(EXPERIMENT,CRYSTAL,MODEL)
% EXPERIMENT contains the scattering geometry
% CRYSTAL - unit cell information
% MODEL - geometry to simulate requires the following:
%   MODEL.UVW_GS
%   MODEL.Uiso_GS
% and if MODEL.useExcitedState
%   MODEL.UVW_ES
%   MODEL.Uiso_ES
% Added in v2:
%   MODEL.fAtomicType = Kirkland OR EDICO
%   MODEL.fAtomicComplex = 1 or 0
%
% NOTE: Vectors are row vectors, so to convert bases use V2 = V1*T12

imgsize = EXPERIMENT.imgsize;
cen = EXPERIMENT.cen;
pixelscale = EXPERIMENT.pixelSize/(EXPERIMENT.detectorDist*EXPERIMENT.magnification);

modk0 = 1/EXPERIMENT.wavelength;

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
HKL_minvals = min(HKL_IMG_L)
HKL_maxvals = max(HKL_IMG_L)
HKL_maxList = createHKL_list(HKL_minvals,HKL_maxvals);
qHKL_maxList = HKL_maxList*CRYSTAL.T_hkl2xyz;
qHKLs_dotk0 = qHKL_maxList*k0_normT;
qIMG_dotk0 = qIMG_L*k0_normT;
qdotk0_min = min(qIMG_dotk0)-EXPERIMENT.qHKLtolerance
qdotk0_max = max(qIMG_dotk0)+EXPERIMENT.qHKLtolerance
indices2keep = (qHKLs_dotk0 >= qdotk0_min) & (qHKLs_dotk0 <= qdotk0_max);
HKL_list = HKL_maxList(indices2keep,:);
qHKL_list = qHKL_maxList(indices2keep,:);
modQ_hkl = sqrt(sum(qHKL_list.^2,2));


%-------------------------------------------------------------------------
% NOTE
% It might be quicker to first remove HKL reflections that don't contribute
% as below, especially if this code is used for fitting. However, tests 
% were quicker without this step
%---------------------------------------------------------------------------
% % Check number of HKL contributing
% nSigma = 2;
% dQtol = nSigma*sqrt((max(modqIMG_L)*MODEL.thetaParallel)^2 + MODEL.sigmaIsotropic^2)
% % 
% % disp('looping HKL');
% % tic
% % [indicesContributing, dQ_nearEwald, mod_dQ_nearEwald] = findContributingReflections(qHKL_list,qIMG_L,dQtol);
% % disp(['Num indicies found = ' num2str(sum(indicesContributing))]);
% % toc
% 
% disp('meshed');
% tic
% [indicesContributing, dQ_nearEwald, mod_dQ_nearEwald] = findContributingReflectionsMeshed(qHKL_list,qIMG_L,dQtol);
% disp(['Num indicies found = ' num2str(sum(indicesContributing))]);
% toc
% 

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

if MODEL.xExcited
    Fhkl_GS = calcFhkl_List(Fatomic, HKL_list, MODEL.UVW_GS, MODEL.Uiso_GS);
    Fhkl_ES = calcFhkl_List(Fatomic, HKL_list, MODEL.UVW_ES, MODEL.Uiso_ES);
    Fhkls = MODEL.xExcited*Fhkl_ES + (1-MODEL.xExcited)*Fhkl_GS;
else
    Fhkls = calcFhkl_List(Fatomic,HKL_list,MODEL.UVW_GS,MODEL.Uiso_GS);
end
Ihkls = abs(Fhkls).^2;
%-------------------------------
% Calculate Diffraction pattern
k0_mesh = meshgrid(k0,1:nPix);
% Basis vectors for decomposing dQ
qBasis_rad = normalizeVectorList(qIMG_L);
qBasis_par = normalizeVectorList(k0_mesh+qIMG_L);
qBasis_rot = cross(qBasis_par,qBasis_rad,2);
% sigma values for describing peak shapes
sigmaRad = MODEL.sigmaIsotropic;
invSigmaRad2rd = 0.5/(sigmaRad^2);
sigmaPar = MODEL.thetaParallel*modqIMG_L+sigmaRad;
invSigmaPar2rd = 0.5./(sigmaPar.^2 + sigmaRad^2);
sigmaRot = MODEL.thetaRotation*modqIMG_L+sigmaRad;
invSigmaRot2rd = 0.5./(sigmaRot.^2 + sigmaRad^2);

nHKL = size(HKL_list,1);
disp(['nHKL = ',num2str(nHKL)]);
IMG_allContrib = zeros([nHKL,nPix]);
Vpixels = 1:nPix;
Normalization = 1./(sigmaRad*sigmaPar.*sigmaRot);

tic
for indx = 1:nHKL
    qHKL_mesh = meshgrid(qHKL_list(indx,:),Vpixels);
    dQ = qHKL_mesh - qIMG_L;
    dQrad2rd = dot(dQ,qBasis_rad,2).^2;
    dQpar2rd = dot(dQ,qBasis_par,2).^2;
    dQrot2rd = dot(dQ,qBasis_rot,2).^2;
    argExp = dQrad2rd*invSigmaRad2rd + dQpar2rd.*invSigmaPar2rd + dQrot2rd.*invSigmaRot2rd;
    IMG_allContrib(indx,:) = Ihkls(indx)*Normalization.*exp(-argExp);
end
toc

IMG_1d = sum(IMG_allContrib,1);
%-----------------------------------------------
% Output
IMG = reshape(IMG_1d,[imgsize.y,imgsize.x]);
if EXPERIMENT.reflectImage
    IMG = fliplr(IMG);
end

%--------------------------------------------------------
%            Nested Functions
%--------------------------------------------------------
    function Fhkl = calcFhkl_List(F_at,HKLs,UVW,Uiso)    
        [Umesh,modQmesh] = meshgrid(Uiso,modQ_hkl);
        FatT = F_at';
        Fhkl_uvw = FatT.*exp(-2*pi*(Umesh.*modQmesh).^2).*cos(2*pi*HKLs*UVW');
        Fhkl = sum(Fhkl_uvw,2);
    end

    function VLnorm = normalizeVectorList(VL)
        modVL = sqrt(sum(VL.^2,2));
        modVLmesh = meshgrid(modVL,1:3)';
        VLnorm = VL./modVLmesh;
    end

    function [hklIndices2use, dQ_contributing, dQ_mod]  = findContributingReflections(qHKL_list,qPix_list,qTolerance)
        % looks for HKLs within a threshold distance from pixels
        % might be quicker to calculate the distance via k0.qHKL but dQ_mod
        % could also be saved for diffraction calc  
        nhkl = size(qHKL_list,1);
        npix = size(qPix_list,1);
        V = 1:npix;
        dQ_all = zeros(nhkl,npix,3);
        for iHKL = 1:nhkl
            qHKL = meshgrid(qHKL_list(iHKL,:),V);
            dQ_all(iHKL,:,:) = qHKL - qPix_list;
        end
        mod_dQ_all = sqrt(sum(dQ_all.^2,3));
        min_mod_dQ_all = min(mod_dQ_all,[],2);
        hklIndices2use = min_mod_dQ_all <= qTolerance;
        dQ_contributing = dQ_all(hklIndices2use,:,:);
        dQ_mod = mod_dQ_all(hklIndices2use,:);
    end


    function [hklIndices2use, dQ_contributing, dQ_mod]  = findContributingReflectionsMeshed(qHKL_list,qPIX_list,qTolerance)
        % looks for HKLs within a threshold distance from pixels
        % might be quicker to calculate the distance via k0.qHKL but dQ_mod
        % is also used for diffraction calc  
        nhkl = size(qHKL_list,1);
        npix = size(qPIX_list,1);
        V = 1:npix;
        qHKL_3d = zeros(nhkl,npix,3);
        qPIX_3d = zeros(nhkl,npix,3);
%         dQ_all = zeros(nhkl,npix,3);

        for iq = 1:3
            [qPIX_3d(:,:,iq),qHKL_3d(:,:,iq)] = meshgrid(qPIX_list(:,iq),qHKL_list(:,iq));
        end
        dQ_all = qHKL_3d - qPIX_3d;
        mod_dQ_all = sqrt(sum(dQ_all.^2,3));
        min_mod_dQ_all = min(mod_dQ_all,[],2);
        hklIndices2use = min_mod_dQ_all <= qTolerance;
        dQ_contributing = dQ_all(hklIndices2use,:,:);
        dQ_mod = mod_dQ_all(hklIndices2use,:);
    end

    function HKL_list = createHKL_list(minhkls,maxhkls)
        nhkl = maxhkls-minhkls+1;
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

    function Tmat = transformMatFrom_k0(k0)
        V3 = k0/norm(k0); % e-beam in z-direction
        Vy = cross(V3,[1,0,0]);
        V2 = Vy/norm(Vy);
        V1 = cross(V2,V3);
        Tmat = [V1;V2;V3];
    end
end