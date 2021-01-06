function refinedStruct = refineStructure_v9(FED, CRYSTAL, MODEL)
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

UVWgs = CRYSTAL.UVW;
% Transformation matrix fractionals -> Cartesian
T_FracCart = CRYSTAL.T_FracCart;
T_CartFrac = inv(T_FracCart);
k0_xyz = MODEL.k0_uvw*T_FracCart;
k0_unit = k0_xyz/norm(k0_xyz);

XYZgs = UVWgs*T_FracCart;
Ugs = CRYSTAL.Uiso;
xExcited = MODEL.xExcited;
HKL = FED.HKL;
HKLt = permute(HKL,[2,1,3]);  
modk = fHKL2modK(HKL,CRYSTAL.axes);
modk2rd = modk.^2;
nSpots = FED.nSpots;
nCont = FED.nConv;
nTimeDelays = length(FED.timeDelays);

%----- Atomic Scattering Factors  --------------------
% Formatted as [atomIndex,contributionIndex,spotIndex]
Fatomic_all = MODEL.Fatomic;
freeAtoms = MODEL.freeAtoms; %Logical
UVWgs_free = UVWgs(freeAtoms,:);
XYZgs_free = XYZgs(freeAtoms,:);
fixedAtoms = ~MODEL.freeAtoms;

%----- Useful variables ---  
nAt_free = sum(MODEL.freeAtoms);
nAt_fixed = CRYSTAL.nAt - nAt_free;
% Renumbering atms -> atms_free
atmVec = 1:CRYSTAL.nAt;
freeAtomNumbers = atmVec(freeAtoms);
atmRenumbFree = 0*(1:CRYSTAL.nAt);
atmRenumbFree(freeAtoms) = 1:nAt_free;
nU_free = sum(MODEL.freeUiso);
freeUiso = MODEL.freeUiso;
fixedUiso = ~freeUiso;
parIndicesFreeX = 1:3*nAt_free;
parIndicesFreeU = (1:nU_free) + parIndicesFreeX(end);
parIndexDW_es = parIndicesFreeX(end)+nU_free+1;
% parIndexDW_gsHot = parIndexDW_es + 1;
Fhkl_0 = zeros(nCont,nSpots);

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
atomsRecalc = freeAtoms | freeUiso;
atomsStatic = ~atomsRecalc;
freeAtoms_recalc = freeAtoms(atomsRecalc);
fixedAtoms_recalc = fixedAtoms(atomsRecalc);
freeUiso_recalc = freeUiso(atomsRecalc);
fixedUiso_recalc = fixedUiso(atomsRecalc);

%----- Atomic scattering factors - subsets ---------
Fatomic_freeAtoms = Fatomic_all(freeAtoms,:,:);
Fatomic_fixedAtoms = Fatomic_all(fixedAtoms,:,:);
Fatomic_recalc = Fatomic_all(atomsRecalc,:,:);
Fatomic_static = Fatomic_all(atomsStatic,:,:);

%----- Atomic DW-factors ------------
DWatomic_all = calcAtomicDWF(Ugs);
DWatomic_static = DWatomic_all(atomsStatic,:,:);
DWatomic_recalc = DWatomic_all(atomsRecalc,:,:);

%----- Atom phases ----
Patomic_all = calcAtomicPhases(UVWgs);
Patomic_recalc = Patomic_all(atomsRecalc,:,:);

%------ GS structure factors --------
Fhkl_GS = calcFhkl(Fatomic_all, UVWgs, Ugs);
Fhkl_static = calcFhkl(Fatomic_static, UVWgs(fixedAtoms,:), Ugs(fixedAtoms));

Igs = Fhkl2Icalc(Fhkl_GS);


%------ Restraint Functions -----------
% calcR = @(X) sqrt(sum(X.^2,2))';

atmPairsRest = MODEL.restraints.atmPairs;
atmPairsRestRenum = atmRenumbFree(atmPairsRest);
DistListGS = atmPairs2dists(atmPairsRest,XYZgs);
if MODEL.restraints.atomPositions
    lambdaAtomCoords = meshgrid(MODEL.restraints.atomPositions,1:3);
    lambdaAtomCoords = reshape(lambdaAtomCoords,[1,3*nAt_free]);
else
    lambdaAtomCoords = Model.restraints.lambdaPos*ones([1,3*nAt_free]);
end

%----- Generate random starting geometries
dXYZ_randomStart = random('norm',0,MODEL.sigmaXYZstart_allT,[nAt_free,3,nTimeDelays]);
dXYZ_randomStart_T1 = random('norm',0,MODEL.sigmaXYZstart_T1,[nAt_free,3]);
dXYZ_randomStart(:,:,1) = dXYZ_randomStart(:,:,1) + dXYZ_randomStart_T1;
DWstart = [0,0];

% Initialize variables for storing results
XYZfree_T = zeros(nAt_free,3,nTimeDelays);
XYZfixed_T = zeros(nAt_fixed,3,nTimeDelays);
XYZ_T = zeros(CRYSTAL.nAt,3,nTimeDelays);
UVW_T = zeros(CRYSTAL.nAt,3,nTimeDelays);
Ufree_T = zeros(nU_free,nTimeDelays);
Uiso_T = zeros(CRYSTAL.nAt,nTimeDelays);

% ------- Start Fitting -----------

for T = 1:nTimeDelays
    fun2minimize = @(P) [(parList2intensityRatios(P)-FED.intensity_ratios(T,:))./FED.meanUncertainties, parList2restraints(P)];
    if MODEL.readXYZstart
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
        options = optimset('Algorithm','levenberg-marquardt','Display','iter','MaxFunEvals',1e5,'TolX',1e-6,'MaxIter',1000);
        if MODEL.useBounds
            nPar = length(Pstart);
            LB = -ones(1,nPar);
            LB(parIndicesFreeAsym) = 0;
            [Pfit, ~, ~, exitflag,output] = lsqnonlin(@(Pfit)fun2minimize(Pfit),Pstart,LB,[],options);
        else
            [Pfit, ~, ~, exitflag,output] = lsqnonlin(@(Pfit)fun2minimize(Pfit),Pstart,[],[],options);
        end
    elseif MODEL.opt_fminunc
        options = optimset('Algorithm','quasi-newton','Display','iter','MaxFunEvals',1e5,'TolX',1e-6,'MaxIter',1000);
%         options = optimset('Algorithm','quasi-newton','HessUpdate','steepdesc','Display','iter','MaxFunEvals',1e5,'TolX',1e-6,'MaxIter',1000);
        [Pfit, ~, exitflag, output] = fminunc(@(Pfit)sum(fun2minimize(Pfit).^2),Pstart,options);
    end
    
    XYZ_T(:,:,T) = XYZgs;
    dXYZfit = parList2positions(Pfit);
    XYZfree_T(:,:,T) = dXYZfit + XYZgs_free;
    XYZ_T(freeAtoms,:,T) = dXYZfit + XYZgs_free;
    XYZfixed_T(:,:,T) = XYZgs(fixedAtoms,:);
    Ufree_T(:,T) = Pfit(parIndicesFreeU);
    Uiso_T(:,T) = Ugs;
    Uiso_T(freeUiso,T) = Ufree_T(:,T);
    RatiosCalc(T,:) = parList2intensityRatios(Pfit);
    Residuals(T,:) = (FED.intensity_ratios(T,:)-parList2intensityRatios(Pfit))./FED.meanUncertainties;
    RestraintVals(T,:) = parList2restraints(Pfit).^2;
    PositionChanges(:,T) = calcR(dXYZfit);
    UVW_T(:,:,T) = XYZ_T(:,:,T)*CRYSTAL.T_CartFrac;
    refinedStruct.dXYZstart(:,:,T) = dXYZstart;
    DW_es(T) = Pfit(parIndexDW_es);
%     DW_gsHot(T) = Pfit(parIndexDW_gsHot);
    Flags(T) = exitflag;
    optimizerOutput{T} = output;
end

XYZ_T(freeAtoms,:,:) = XYZfree_T;
XYZ_T(fixedAtoms,:,:) = XYZfixed_T;

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
% refinedStruct.DW_gsHot = DW_gsHot;
refinedStruct.exitFlags = Flags;


    function DW_at = calcAtomicDWF(U)
        DW_at = zeros([length(U),size(modk2rd,1),size(modk2rd,2)]);
        for S = 1:nSpots
            [modk2rdmesh, Umesh] = meshgrid(modk2rd(:,S),U);
            DW_at(:,:,S) = exp(-2*pi^2*Umesh.*modk2rdmesh);
        end
    end

    function P_at = calcAtomicPhases(uvw)
        P_at = zeros([size(uvw,1),nCont,nSpots]);
        for S = 1:nSpots
            P_at(:,:,S) = cos(2*pi*uvw*HKLt(:,:,S));
        end
    end

% Old version (v4) - Encountered a problem in v2017a sum() when nConv = 1
%     function Fhkl = calcFhkl(F_at,uvw,U)
%         Fhkl = zeros(nConv,nSpots);       
%         for S = 1:nSpots
%             [modkmesh, Umesh] = meshgrid(modk(:,S),U);
%             F_uvw_hkl = F_at(:,:,S).*exp(-2*pi*(Umesh.*modkmesh).^2).*cos(2*pi*uvw*HKLt(:,:,S));
%             Fhkl(:,S) = sum(F_uvw_hkl,1);
%         end
%     end

    function Fhkl = calcFhkl(F_at,uvw,U)
        Fhkl = zeros(nCont,nSpots);       
        for S = 1:nSpots
            for C = 1:nCont
                F_uvw_hkl = F_at(:,C,S).*exp(-2*pi^2*U*modk2rd(C,S)).*cos(2*pi*uvw*HKLt(:,C,S));
                Fhkl(C,S) = sum(F_uvw_hkl);
            end
        end
    end

    function Fhkl = calcAllFhkl(uvw,U)
        Fhkl = zeros(nCont,nSpots);       
        for S = 1:nSpots
            for C = 1:nCont
                F_uvw_hkl = Fatomic_all(:,C,S).*exp(-2*pi^2*U*modk2rd(C,S)).*cos(2*pi*uvw*HKLt(:,C,S));
                Fhkl(C,S) = sum(F_uvw_hkl);
            end
        end
    end

    function Fhkl = recalcFhkl(uvw,U)
        Fhkl = zeros(nCont,nSpots);
        for S = 1:nSpots
            [modk2rdmesh, Umesh] = meshgrid(modk2rd(:,S).^2,U);
            B = DWatomic_recalc(:,:,S);
            B(freeUiso_recalc,:) = exp(-2*pi^2*Umesh.*modk2rdmesh);
            P = Patomic_recalc(:,:,S);
            P(freeAtoms_recalc,:) = cos(2*pi*uvw*HKLt(:,:,S));
            F = Fatomic_recalc(:,:,S).*B.*P;
            Fhkl(:,S) = sum(F,1);
        end
    end

    function Icalc = Fhkl2Icalc(Fcalc)
        Icalc = sum(FED.Wts.*Fcalc.^2,1);
    end

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

    function Ratios = parList2intensityRatios(PAR)
        UVWes = UVWgs;
        Ues = Ugs;
        dXYZ = reshape(PAR(parIndicesFreeX)*T_par2dx,[3,nAt_free])';
        UVWes(freeAtoms,:) = dXYZ*T_CartFrac + UVWgs_free;
        Ues(freeUiso) = PAR(parIndicesFreeU);
        Ues = Ues + PAR(parIndexDW_es);
        Fhkl_ES = calcAllFhkl(UVWes,Ues);
        Fhkl_FED = xExcited*Fhkl_ES + (1-xExcited)*Fhkl_GS;
        Ratios = Fhkl2Icalc(Fhkl_FED)./Igs;
    end

    function Restraints = parList2restraints(PAR)
        dXYZ = parList2positions(PAR);
        dXdotK = k0_unit*dXYZ';
        XYZ = dXYZ + XYZgs_free;
        dR = (atmPairs2dists(atmPairsRestRenum,XYZ) - DistListGS)';
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