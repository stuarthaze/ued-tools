function refinedStruct = refineStructure_v6_fminunc(FED, CRYSTAL, MODEL)
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
nSpots = FED.nSpots;
nCont = FED.nConv;
nTimeDelays = length(FED.timeDelays);

%----- Atomic Scattering Factors  --------------------
% Formatted as [atomIndex,contributionIndex,spotIndex]
Fatomic_all = MODEL.Fatomic;
freeAtoms = MODEL.freeAtoms;
UVWgs_free = UVWgs(freeAtoms,:);
XYZgs_free = XYZgs(freeAtoms,:);
fixedAtoms = ~MODEL.freeAtoms;

%----- Useful variables ---  
nAt_free = sum(MODEL.freeAtoms);
nAt_fixed = CRYSTAL.nAt - nAt_free;  
nU_free = sum(MODEL.freeUiso);
freeUiso = MODEL.freeUiso;
fixedUiso = ~freeUiso;
parIndicesFreeX = 1:3*nAt_free;
parIndicesFreeU = (1:nU_free) + parIndicesFreeX(end);
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
T_dx2par = inv(T_par2dx);

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

%----- Atomic B-factors ------------
Batomic_all = calcAtomicBfactors(Ugs);
Batomic_static = Batomic_all(atomsStatic,:,:);
Batomic_recalc = Batomic_all(atomsRecalc,:,:);

%----- Atom phases ----
Patomic_all = calcAtomicPhases(UVWgs);
Patomic_recalc = Patomic_all(atomsRecalc,:,:);

%------ GS structure factors --------
Fhkl_GS = calcFhkl(Fatomic_all, UVWgs, Ugs);
Fhkl_static = calcFhkl(Fatomic_static, UVWgs(fixedAtoms,:), Ugs(fixedAtoms));

Igs = Fhkl2Icalc(Fhkl_GS);


%------ Restraint Functions -----------
fCalc_dR = @(X) sqrt(sum(X.^2,2))';
% A = 0.02; % Distance for transition between quadratic and linear penalty
% fPositionPenaltyLinear = @(X) abs(sqrt(exp(-(X/A).^2).*(X.^2)/(2*A)+(1-exp(-(X/A).^2)).*(abs(X)-0.5*A)));
% fRestraintFunction = @(X) MODEL.restraints.lambdaPos*fPositionPenaltyLinear(fCalc_dR(fPositionChanges(parList2positions(X))));
if MODEL.restraints.lambdaPos && ~(MODEL.restraints.lambdaBeamDirection)
    fRestraintFunction = @(PX) MODEL.restraints.lambdaPos*PX(parIndicesFreeX)*T_par2dx;
elseif ~MODEL.restraints.lambdaPos && MODEL.restraints.lambdaBeamDirection
    fRestraintFunction = @(PX) MODEL.restraints.lambdaBeamDirection*dxInBeamDirection(PX);
elseif MODEL.restraints.lambdaPos && MODEL.restraints.lambdaBeamDirection
    fRestraintFunction = @(PX) [MODEL.restraints.lambdaPos*PX(parIndicesFreeX)*T_par2dx, MODEL.restraints.lambdaBeamDirection*dxInBeamDirection(PX)];
else
    fRestraintFunction = @(PX) [];
end
    

%----- Generate random starting geometries
dXYZ_randomStart = random('norm',0,MODEL.sigmaXYZstart_allT,[nAt_free,3,nTimeDelays]);
dXYZ_randomStart_T1 = random('norm',0,MODEL.sigmaXYZstart_T1,[nAt_free,3]);
dXYZ_randomStart(:,:,1) = dXYZ_randomStart(:,:,1) + dXYZ_randomStart_T1;

% Initialize variables for storing results
XYZfree_T = zeros(nAt_free,3,nTimeDelays);
XYZfixed_T = zeros(nAt_fixed,3,nTimeDelays);
XYZ_T = zeros(CRYSTAL.nAt,3,nTimeDelays);
UVW_T = zeros(CRYSTAL.nAt,3,nTimeDelays);
Ufree_T = zeros(nU_free,nTimeDelays);
Uiso_T = zeros(CRYSTAL.nAt,nTimeDelays);

% ------- Start Fitting -----------

for T = 1:nTimeDelays
    fun2minimize = @(P) [(parList2intensityRatios(P)-FED.intensity_ratios(T,:))./FED.meanUncertainties, fRestraintFunction(P)];
    if T == 1
        dXYZstart = dXYZ_randomStart(:,:,1);        
    elseif MODEL.useSequentialStart
        dXYZstart = dXYZ_randomStart(:,:,T) + dXYZfit;
    else
        dXYZstart = dXYZ_randomStart(:,:,T);
    end
    Pstart = [reshape(dXYZstart',[1,3*nAt_free])*T_dx2par,Ugs(MODEL.freeUiso)'];
    
    if ~MODEL.opt_fminunc
        options = optimset('Algorithm','levenberg-marquardt','Display','iter','MaxFunEvals',1e5,'TolX',1e-6,'MaxIter',500);
        if MODEL.useBounds
            nPar = length(Pstart);
            LB = -ones(1,nPar);
            LB(parIndicesFreeAsym) = 0;
            Pfit = lsqnonlin(@(Pfit)fun2minimize(Pfit),Pstart,LB,[],options);
        else
            Pfit = lsqnonlin(@(Pfit)fun2minimize(Pfit),Pstart,[],[],options);
        end
    elseif MODEL.opt_fminunc
        options = optimset('Algorithm','quasi-newton','Display','iter','MaxFunEvals',1e5,'TolX',1e-6,'MaxIter',500);
        Pfit = fminunc(@(Pfit)sum(fun2minimize(Pfit).^2),Pstart,options);
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
    RestraintVals(T,:) = fRestraintFunction(Pfit).^2;
    PositionChanges(:,T) = fCalc_dR(dXYZfit);
    UVW_T(:,:,T) = XYZ_T(:,:,T)*CRYSTAL.T_CartFrac;
    refinedStruct.dXYZstart(:,:,T) = dXYZstart;
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


    function B_at = calcAtomicBfactors(U)
        B_at = zeros([length(U),size(modk,1),size(modk,2)]);
        for S = 1:nSpots
            [modkmesh, Umesh] = meshgrid(modk(:,S),U);
            B_at(:,:,S) = exp(-2*pi*(Umesh.*modkmesh).^2);
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
                F_uvw_hkl = F_at(:,C,S).*exp(-2*pi*(U*modk(C,S)).^2).*cos(2*pi*uvw*HKLt(:,C,S));
                Fhkl(C,S) = sum(F_uvw_hkl);
            end
        end
    end

    function Fhkl = recalcFhkl(uvw,U)
        Fhkl = zeros(nCont,nSpots);
        for S = 1:nSpots
            [modkmesh, Umesh] = meshgrid(modk(:,S),U);
            B = Batomic_recalc(:,:,S);
            B(freeUiso_recalc,:) = exp(-2*pi*(Umesh.*modkmesh).^2);
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

    function Ratios = parList2intensityRatios(PAR)
        dXYZ = reshape(PAR(parIndicesFreeX)*T_par2dx,[3,nAt_free])';
        UVWfree = dXYZ*T_CartFrac + UVWgs_free;
        Ufree = PAR(parIndicesFreeU);
        Fhkl_free = recalcFhkl(UVWfree,Ufree);
        Fhkl_FED = xExcited*(Fhkl_static + Fhkl_free) + (1-xExcited)*Fhkl_GS;
        Ratios = Fhkl2Icalc(Fhkl_FED)./Igs;
    end

    function X = parList2positions(PAR)
        X = reshape(PAR(parIndicesFreeX)*T_par2dx,[3,nAt_free])';
    end
    
    function dX = dxInBeamDirection(PAR)
        dX = k0_unit*reshape(PAR(parIndicesFreeX)*T_par2dx,[3,nAt_free]);
    end
end