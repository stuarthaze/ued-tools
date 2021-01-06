function refinedStruct = refineStructure_v13(FED, CRYSTAL, MODEL)
% Changes: 
% v3
% Use changes in position (dxyz) rather than coordinates (xyz) as parameters
% Restrain motion along beam direction
% v4
% Use random displacements from GS as starting points
% v5
% Moved to Matlab 2017a and needed to change calcFhkl()
% v6
% Option for different coordinate systems.
% v7
% Addition of bond-length restraints. Corrected DW-factor. 
% v8
% Added global DWF
% v9
% Added individual atomic position restraints: MODEL.restraints.atomPositions
% v10
% New way to calculate Iroi = Wts*Ihkl
% v11
% Sped up version
% v12
% Added option for using DWF
% v13
% Option to restrain atom positions to specified positions

UVWgs = CRYSTAL.UVW;
% UVWes = UVWgs;
% Transformation matrix fractionals -> Cartesian
T_FracCart = CRYSTAL.T_FracCart;
T_CartFrac = inv(T_FracCart);
k0_xyz = MODEL.k0_uvw*T_FracCart;
k0_unit = k0_xyz/norm(k0_xyz);

XYZgs = UVWgs*T_FracCart;
XYZtemp = XYZgs;
Ugs = CRYSTAL.Uiso;
xExcited = MODEL.xExcited;
HKL = FED.HKL;
% nHKL = size(HKL,1);
qHKL = HKL*CRYSTAL.T_hkl2xyz;
modQhkl2rd = sum(qHKL.^2,2)';
WTS = FED.Wts;
HKLt = HKL';
% modk = fHKL2modK(HKL,CRYSTAL.axes);
% % modk2rd = modk.^2;
% nROIs = FED.nROIs; 
% nConv = FED.nConv;
nTimeDelays = length(FED.timeDelays);

%----- Atomic Scattering Factors  --------------------
% Formatted as [atomIndex,HKL_index]

% Calculate Fatomic for HKLs
if strcmp(MODEL.fAtomicType,'Kirkland')
    Fatomic = calculateAtomicF_Kirkland(CRYSTAL,FED.HKL);
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

% Fatomic = MODEL.Fatomic;
freeAtoms = MODEL.freeAtoms; %Logical
UVWgs_free = UVWgs(freeAtoms,:);
XYZgs_free = XYZgs(freeAtoms,:);
fixedAtoms = ~MODEL.freeAtoms;

%----- Useful variables ---  
nAt_free = sum(MODEL.freeAtoms);
nAt_fixed = CRYSTAL.nAt - nAt_free;
% Renumbering atms -> atms_free
atmVec = 1:CRYSTAL.nAt;
% freeAtomNumbers = atmVec(freeAtoms);
% atmRenumbFree = 0*(1:CRYSTAL.nAt);
% atmRenumbFree(freeAtoms) = 1:nAt_free;
nU_free = sum(MODEL.freeUiso);
freeUiso = MODEL.freeUiso;
% fixedUiso = ~freeUiso;
parIndicesFreeX = 1:3*nAt_free;
parIndicesFreeU = (1:nU_free) + 3*nAt_free;

if MODEL.refineDWF
    parIndexDW_es = 3*nAt_free + nU_free+1;
    DWstart = 0;
else
    parIndexDW_es = [];
    DWstart = [];
end
% parIndexDW_gsHot = parIndexDW_es + 1;
% Fhkl_0 = zeros(nConv,nROIs);

% Transformation matrix Parameters -> dX_vector
if MODEL.useSymmetrizedCoords == 1
    nX_free_asym = 3*nAt_free/MODEL.nSymCoords;
    T_par2dx = [eye(nX_free_asym), eye(nX_free_asym); ...
                eye(nX_free_asym),-eye(nX_free_asym)];
    parIndicesFreeAsym = 1:nX_free_asym;
else
    T_par2dx = eye(3*nAt_free);
end
T_dx2par = eye(3*nAt_free)/T_par2dx;

%----- Atom indicies for recalculating ---
% atomsRecalc = freeAtoms | freeUiso;
% atomsStatic = ~atomsRecalc;
% freeAtoms_recalc = freeAtoms(atomsRecalc);
% fixedAtoms_recalc = fixedAtoms(atomsRecalc);
% freeUiso_recalc = freeUiso(atomsRecalc);
% fixedUiso_recalc = fixedUiso(atomsRecalc);

% %----- Atomic scattering factors - subsets ---------
% Fatomic_freeAtoms = Fatomic(freeAtoms,:);
% Fatomic_fixedAtoms = Fatomic(fixedAtoms,:);
% Fatomic_recalc = Fatomic(atomsRecalc,:);
% Fatomic_static = Fatomic(atomsStatic,:);

%----- Atomic DW-factors ------------
DW_at_GS = calcAtomicDWF(Ugs);
% DW_at = DW_at_GS;
% DWatomic_static = DWatomic_all(atomsStatic,:);
% DWatomic_recalc = DWatomic_all(atomsRecalc,:);

%----- Atom phases ----
P_at_GS = calcAtomicPhases(UVWgs);
% P_at = P_at_GS;
% Patomic_recalc = Pat_GS(atomsRecalc,:);

%------ GS structure factors --------
F_uvwhkl_GS = Fatomic.*DW_at_GS.*P_at_GS;
% F_uvwhkl = F_uvwhkl_GS;
Fhkl_GS_new = sum(F_uvwhkl_GS,1);
Fhkl_GS = calcFhkl(Fatomic, UVWgs, Ugs);
maxFhklDiff = max(Fhkl_GS - Fhkl_GS_new);
% Fhkl_static = calcFhkl(Fatomic_static, UVWgs(fixedAtoms,:), Ugs(fixedAtoms));

Igs = Fhkl2Icalc(Fhkl_GS);


%------ Restraint Functions -----------
% calcR = @(X) sqrt(sum(X.^2,2))';

% atmPairsRestRenum = atmRenumbFree(atmPairsRest);
DistListGS = atmPairs2dists(MODEL.restraints.atomPairs,XYZgs);
lambdaAtomPositionsFree = MODEL.restraints.atomPositions(MODEL.freeAtoms);
lambdaAtomCoords3xN = meshgrid(lambdaAtomPositionsFree,1:3);
lambdaAtomCoords = reshape(lambdaAtomCoords3xN,[1,3*nAt_free]);

%----- Generate random starting geometries
dXYZ_randomStart = random('norm',0,MODEL.sigmaXYZstart_allT,[nAt_free,3,nTimeDelays]);
dXYZ_randomStart_T1 = random('norm',0,MODEL.sigmaXYZstart_T1,[nAt_free,3]);
dXYZ_randomStart(:,:,1) = dXYZ_randomStart(:,:,1) + dXYZ_randomStart_T1;
% DWstart = [0,0];

% Initialize variables for storing results
XYZfree_T = zeros(nAt_free,3,nTimeDelays);
XYZfixed_T = zeros(nAt_fixed,3,nTimeDelays);
XYZ_T = zeros(CRYSTAL.nAt,3,nTimeDelays);
dXYZ_T = zeros(nAt_free,3,nTimeDelays);
UVW_T = zeros(CRYSTAL.nAt,3,nTimeDelays);
Ufree_T = zeros(nU_free,nTimeDelays);
Uiso_T = zeros(CRYSTAL.nAt,nTimeDelays);

% ------- Start Fitting -----------

for T = 1:nTimeDelays
    if MODEL.restraints.useXYZT
        XYZrest = MODEL.restraints.XYZT(freeAtoms,:,T);
        fun2minimize = @(P) [(parList2intensityRatios(P)-FED.ratios(T,:))./FED.meanUncertainties_R, parList2restraints_nonZeroXYZ(P,XYZrest)];
    else
        fun2minimize = @(P) [(parList2intensityRatios(P)-FED.ratios(T,:))./FED.meanUncertainties_R, parList2restraints(P)];
    end
    if MODEL.readXYZTstart
        if size(MODEL.XYZTstart,3) ~= FED.nTimeDelays
            disp('Error: number of start points different from number of time points');
            return
        end
        dXYZstart = MODEL.XYZTstart(freeAtoms,:,T) - XYZgs_free;
    elseif MODEL.readXYZstart
        dXYZstart = MODEL.XYZstart(freeAtoms,:) - XYZgs_free;
    elseif T == 1
        dXYZstart = dXYZ_randomStart(:,:,1);
    elseif MODEL.useSequentialStart
        dXYZstart = dXYZ_randomStart(:,:,T) + dXYZfit;
    else
        dXYZstart = dXYZ_randomStart(:,:,T);
    end
    Pstart = [reshape(dXYZstart',[1,3*nAt_free])*T_dx2par,Ugs(MODEL.freeUiso)',DWstart];
    
    if ~MODEL.opt_fminunc
        options = optimset('Algorithm','levenberg-marquardt','Display','iter','MaxFunEvals',5e5,'TolX',1e-5,'MaxIter',1000);
        if MODEL.useBounds
            nPar = length(Pstart);
            LB = -ones(1,nPar);
            LB(parIndicesFreeAsym) = 0;
            [Pfit, ~, ~, exitflag,output] = lsqnonlin(@(Pfit)fun2minimize(Pfit),Pstart,LB,[],options);
        else
            [Pfit, ~, ~, exitflag,output] = lsqnonlin(@(Pfit)fun2minimize(Pfit),Pstart,[],[],options);
        end
    elseif MODEL.opt_fminunc
        options = optimset('Algorithm','quasi-newton','Display','iter','MaxFunEvals',2e5,'MaxIter',1000);
%         options = optimset('Algorithm','quasi-newton','HessUpdate','steepdesc','Display','iter','MaxFunEvals',1e5,'TolX',1e-6,'MaxIter',1000);
        [Pfit, ~, exitflag, output] = fminunc(@(Pfit)sum(fun2minimize(Pfit).^2),Pstart,options);
    end
    
    XYZ_T(:,:,T) = XYZgs;
    dXYZfit = parList2positions(Pfit);
    dXYZ_T(:,:,T) = dXYZfit;
    XYZfree_T(:,:,T) = dXYZfit + XYZgs_free;
    XYZ_T(freeAtoms,:,T) = dXYZfit + XYZgs_free;
%     XYZfixed_T(:,:,T) = XYZgs(fixedAtoms,:);
    Ufree_T(:,T) = Pfit(parIndicesFreeU);
    Uiso_T(:,T) = Ugs;
    Uiso_T(freeUiso,T) = Ufree_T(:,T);
    RatiosCalc(T,:) = parList2intensityRatios(Pfit);
    Residuals(T,:) = (FED.ratios(T,:)-parList2intensityRatios(Pfit))./FED.meanUncertainties_R;
    RestraintVals(T,:) = parList2restraints(Pfit).^2;
    PositionChanges(:,T) = calcR(dXYZfit);
    UVW_T(:,:,T) = XYZ_T(:,:,T)*CRYSTAL.T_CartFrac;
    refinedStruct.dXYZstart(:,:,T) = dXYZstart;
    if MODEL.refineDWF
        DW_es(T) = Pfit(parIndexDW_es);
    else
        DW_es(T) = 0;
    end
%     DW_gsHot(T) = Pfit(parIndexDW_gsHot);
    Flags(T) = exitflag;
    optimizerOutput{T} = output;
end


refinedStruct.dXYZ_T = dXYZ_T;
refinedStruct.XYZ_T = XYZ_T;
refinedStruct.UVW_T = UVW_T;
refinedStruct.XYZfree_T = XYZfree_T;
refinedStruct.Ufree_T = Ufree_T;
refinedStruct.RatiosCalc = RatiosCalc;
refinedStruct.initialResiduals = FED.initialResiduals;
refinedStruct.Residuals = Residuals;
refinedStruct.RestraintVals = RestraintVals;
refinedStruct.PositionChanges = PositionChanges;
refinedStruct.timeDelays = FED.timeDelays;
refinedStruct.Uiso_T = Uiso_T;
refinedStruct.DW_es = DW_es;
refinedStruct.atomNumbersRefining = atmVec(MODEL.freeAtoms);
% refinedStruct.DW_gsHot = DW_gsHot;
refinedStruct.exitFlags = Flags;


    function DW_at = calcAtomicDWF(U)
        [modQ2rdmesh, Umesh] = meshgrid(modQhkl2rd,U);
        DW_at = exp(-2*pi^2*Umesh.*modQ2rdmesh);
    end

    function Pout = calcAtomicPhases(uvw)
        Pout = cos(2*pi*uvw*HKLt);
    end

    function Fhkl = calcFhkl(F_at,uvw,U)
        [modQ2rdmesh, Umesh] = meshgrid(modQhkl2rd,U);
        DW = exp(-2*pi^2*Umesh.*modQ2rdmesh);
        P = cos(2*pi*uvw*HKLt);
        F_uvw_hkl = F_at.*DW.*P;
        Fhkl = sum(F_uvw_hkl);
    end

    function Icalc = Fhkl2Icalc(Fcalc)
        Icalc = (Fcalc.^2)*WTS;
    end

    function Ratios = parList2intensityRatios(PAR)
        P_at = P_at_GS;
        DW_at = DW_at_GS;
        dXYZ = reshape(PAR(parIndicesFreeX)*T_par2dx,[3,nAt_free])';
        UVWfree =  dXYZ*T_CartFrac + UVWgs_free;
        P_at(freeAtoms,:) = calcAtomicPhases(UVWfree);
        Ues = PAR(parIndicesFreeU);
        DW_at(freeUiso,:) = calcAtomicDWF(Ues);
        F_uvwhkl_ES = Fatomic.*DW_at.*P_at;
        if MODEL.refineDWF
            DW_ES = exp(-2*pi^2*PAR(parIndexDW_es)*modQhkl2rd);
            Fhkl_ES = sum(F_uvwhkl_ES,1).*DW_ES;
        else
            Fhkl_ES = sum(F_uvwhkl_ES,1);
        end
        Fhkl_FED = xExcited*Fhkl_ES + (1-xExcited)*Fhkl_GS;
        Ratios = Fhkl2Icalc(Fhkl_FED)./Igs;
    end

%     function Fhkl = calcAllFhkl(uvw,U)
%         [modQ2rdmesh, Umesh] = meshgrid(modQhkl2rd,U);
%         DW_at = exp(-2*pi^2*Umesh.*modQ2rdmesh);
%         P_at = cos(2*pi*uvw*HKLt);
%         F_uvw_hkl = Fatomic.*DW_at.*P_at;
%         Fhkl = sum(F_uvw_hkl);
%     end
% 
%     function Fhkl = recalcFhkl(uvw,U)
%         [modQhkl2rdmesh, Umesh] = meshgrid(modQhkl2rd.^2,U);
%         B = DWatomic_recalc;
%         B(freeUiso_recalc,:) = exp(-2*pi^2*Umesh.*modQhkl2rdmesh);
%         P = Patomic_recalc;
%         P(freeAtoms_recalc,:) = cos(2*pi*uvw*HKLt(:,:,S));
%         F = Fatomic_recalc.*B.*P;
%         Fhkl(:,S) = sum(F,1);
%     end

%     function Ratios = calcRatios(XYZfree,Ufree)
%         Ratios = zeros(1,nSpots);
%         Fhkl_free = recalcFhkl(XYZfree,Ufree);
%         Fhkl_FED = xExcited*(Fhkl_static + Fhkl_free) + (1-xExcited)*Fhkl_GS;
%         Ratios = Fhkl2Icalc(Fhkl_FED)./Igs;
%     end
%         
%     function Ratios = parList2intensityRatios(PAR)
%         dXYZ = reshape(PAR(parIndicesFreeX),[3,nAt_free])';
%         UVWfree = dXYZ*T_CartFrac + UVWgs_free;
%         Ufree = PAR(parIndicesFreeU);
%         Fhkl_free = recalcFhkl(UVWfree,Ufree);
%         Fhkl_FED = xExcited*(Fhkl_static + Fhkl_free) + (1-xExcited)*Fhkl_GS;
%         Ratios = Fhkl2Icalc(Fhkl_FED)./Igs;
%     end


    function X = parList2positions(PAR)
        X = reshape(PAR(parIndicesFreeX)*T_par2dx,[3,nAt_free])';
    end
    
    function dX = dxInBeamDirection(PAR)
        dX = k0_unit*reshape(PAR(parIndicesFreeX)*T_par2dx,[3,nAt_free]);
    end

    function R = calcR(XYZ)
        R = sqrt(sum(XYZ.^2,2));
    end

    function Rab = atmPairs2dists(AB_list,XYZ)
        A = AB_list(:,1);
        B = AB_list(:,2);
        Rab = calcR(XYZ(A,:)-XYZ(B,:));
    end
% 
%     function Ratios = parList2intensityRatios(PAR)
%         dXYZ = reshape(PAR(parIndicesFreeX)*T_par2dx,[3,nAt_free])';
%         UVWfree = dXYZ*T_CartFrac + UVWgs_free;
%         Ufree = PAR(parIndicesFreeU);
%         DWF_es = exp(-2*pi^2*modk2rd*PAR(parIndexDW_es));
%         DWF_gsHot = exp(-2*pi^2*modk2rd*PAR(parIndexDW_gsHot));
%         Fhkl_free = recalcFhkl(UVWfree,Ufree);
%         Fhkl_FED = xExcited*(Fhkl_static + Fhkl_free).*DWF_es + (1-xExcited)*Fhkl_GS.*DWF_gsHot;
%         Ratios = Fhkl2Icalc(Fhkl_FED)./Igs;
%     end


    function Restraints = parList2restraints_nonZeroXYZ(PAR,XYZrestraint)
        
        dXYZ = parList2positions(PAR);
%         dXdotK = k0_unit*dXYZ'; Use to restrain in beam direction
        XYZall = XYZgs;
        XYZall(freeAtoms,:) = dXYZ + XYZgs_free;
        dX = reshape((dXYZ + XYZgs_free - XYZrestraint)',[1,3*nAt_free]);
        dR = (atmPairs2dists(MODEL.restraints.atomPairs,XYZall) - DistListGS)';
        Restraints = [lambdaAtomCoords.*dX, MODEL.restraints.lambdaDists*dR];
    end


    function Restraints = parList2restraints(PAR)
        dXYZ = parList2positions(PAR);
%         dXdotK = k0_unit*dXYZ'; Use to restrain in beam direction
        XYZall = XYZgs;
        XYZall(freeAtoms,:) = dXYZ + XYZgs_free;
        dR = (atmPairs2dists(MODEL.restraints.atomPairs,XYZall) - DistListGS)';
        Restraints = [lambdaAtomCoords.*PAR(parIndicesFreeX)*T_par2dx, MODEL.restraints.lambdaDists*dR];
    end
% Old version with restraint in beam direction
%     function Restraints = parList2restraints(PAR)
%         dXYZ = parList2positions(PAR);
%         dXdotK = k0_unit*dXYZ';
%         XYZ = dXYZ + XYZgs_free;
%         dR = (atmPairs2dists(atmPairsRestRenum,XYZ) - DistListGS)';
%         Restraints = [lambdaAtomCoords.*PAR(parIndicesFreeX)*T_par2dx,...
%                       MODEL.restraints.lambdaBeamDirection*dXdotK,...
%                       MODEL.restraints.lambdaDists*dR];
%     end

end
% 
% function R = cart2sphPolar(XYZ)
%     % Input: XYZ is an array [nAtom x 3]
%     % OUTPUT:
%     % R.r = radius
%     % R.az = azimuthal is atan(y/x)
%     % R.el = elevation is angle from xy plane = pi/2 - theta
%     [R.az,R.el,R.r] = cart2sph(XYZ(:,1),XYZ(:,2),XYZ(:,3));
% end
% 
% function XYZ = sphPolar2cart(R)
%     [x,y,z] = sph2cart(R.az,R.el,R.r);
%     XYZ = [x,y,z];
% end