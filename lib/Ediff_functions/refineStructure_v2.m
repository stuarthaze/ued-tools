function refinedStruct = refineStructure_v2(FED, CRYSTAL, MODEL)

XYZgs = CRYSTAL.XYZ;
XYZgs_cart = fractionals2cartesians(XYZgs,CRYSTAL.axes);
Ugs = CRYSTAL.Uiso;
xExcited = MODEL.xExcited;
HKL = FED.HKL;
HKLt = permute(HKL,[2,1,3]);
modk = fHKL2modK(HKL,CRYSTAL.axes);
nSpots = FED.nSpots;
nConv = FED.nConv;
nTimeDelays = length(FED.timeDelays);



%----- Atomic Scattering Factors  --------------------
% Formatted as [atomIndex,contributionIndex,spotIndex]
Fatomic_all = MODEL.Fatomic;
freeAtoms = MODEL.freeAtoms;
fixedAtoms = ~MODEL.freeAtoms;

%----- Useful variables ---
nAt_free = sum(MODEL.freeAtoms);
nAt_fixed = CRYSTAL.nAt - nAt_free;
nU_free = sum(MODEL.freeUiso);
freeUiso = MODEL.freeUiso;
fixedUiso = ~freeUiso;
parIndicesXYZ = 1:3*nAt_free;
parIndicesU = (1:nU_free) + parIndicesXYZ(end);
Fhkl_0 = zeros(nConv,nSpots);

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
Patomic_all = calcAtomicPhases(XYZgs);
Patomic_recalc = Patomic_all(atomsRecalc,:,:);

%------ GS structure factors --------
Fhkl_GS = calcFhkl(Fatomic_all, XYZgs, Ugs);
Fhkl_static = calcFhkl(Fatomic_static, XYZgs(fixedAtoms,:), Ugs(fixedAtoms));

Igs = Fhkl2Icalc(Fhkl_GS);

%------ Restraint Functions -----------
fCalc_dR = @(X) sqrt(sum(X.^2,2))';
A = 0.02; % Distance for transition between quadratic and linear penalty
fPositionPenaltyLinear = @(X) abs(sqrt(exp(-(X/A).^2).*(X.^2)/(2*A)+(1-exp(-(X/A).^2)).*(abs(X)-0.5*A)));
fPositionChanges = @(X) fractionals2cartesians(X,CRYSTAL.axes)-XYZgs_cart(MODEL.freeAtoms,:);
fRestraintFunction = @(X) MODEL.restraints.lambdaPos*fPositionPenaltyLinear(fCalc_dR(fPositionChanges(parList2positions(X))));
if MODEL.displayRestraintFunction
    xTest = -0.5:0.01:0.5;
    yTest = fPositionPenaltyLinear(xTest).^2;
    figure(); plot(xTest,yTest);
    title('Position Penalty Function');
end

%------- Start Fitting -----------
Pstart = [reshape(XYZgs(MODEL.freeAtoms,:)',[1,3*nAt_free]),Ugs(MODEL.freeUiso)'];
XYZfree_T = zeros(nAt_free,3,nTimeDelays);
XYZfixed_T = zeros(nAt_fixed,3,nTimeDelays);

Ufree_T = zeros(nU_free,nTimeDelays);
options = optimset('Display','iter','MaxFunEvals',1e5,'TolX',1e-6,'MaxIter',500);
for T = 1:nTimeDelays
    fun2minimize = @(P) [(parList2intensityRatios(P)-FED.intensity_ratios(T,:))./FED.meanUncertainties, fRestraintFunction(P)];
    Pfit = lsqnonlin(@(Pfit)fun2minimize(Pfit),Pstart,[],[],options);
    XYZfit = parList2positions(Pfit);
    XYZfree_T(:,:,T) = XYZfit;
    XYZfixed_T(:,:,T) = XYZgs(fixedAtoms,:);
    Ufree_T(:,T) = Pfit(parIndicesU);
    Residuals(T,:) = (FED.intensity_ratios(T,:)-parList2intensityRatios(Pfit))./FED.meanUncertainties;
    RestraintVals(T,:) = fRestraintFunction(Pfit);
    PositionChanges(:,T) = fCalc_dR(fPositionChanges(XYZfit));
    if MODEL.startFromLastFit
        Pstart = Pfit;
    end
end

XYZ_T(freeAtoms,:,:) = XYZfree_T;
XYZ_T(fixedAtoms,:,:) = XYZfixed_T;

refinedStruct.XYZ_T = XYZ_T;
refinedStruct.XYZfree_T = XYZfree_T;
refinedStruct.XYZfixed_T = XYZfixed_T;
refinedStruct.Ufree_T = Ufree_T;
refinedStruct.Residuals = Residuals;
refinedStruct.RestraintVals = RestraintVals;
refinedStruct.PositionChanges = PositionChanges;
refinedStruct.timeDelays = FED.timeDelays;


    function B_at = calcAtomicBfactors(U)
        B_at = zeros([length(U),size(modk,1),size(modk,2)]);
        for S = 1:nSpots
            [modkmesh, Umesh] = meshgrid(modk(:,S),U);
            B_at(:,:,S) = exp(-2*pi*(Umesh.*modkmesh).^2);
        end
    end

    function P_at = calcAtomicPhases(xyz)
        P_at = zeros([size(xyz,1),nConv,nSpots]);
        for S = 1:nSpots
            P_at(:,:,S) = cos(2*pi*xyz*HKLt(:,:,S));
        end
    end

    function Fhkl = calcFhkl(F_at,xyz,U)
        Fhkl = zeros(nConv,nSpots);       
        for S = 1:nSpots
            [modkmesh, Umesh] = meshgrid(modk(:,S),U);
            F_xyz_hkl = F_at(:,:,S).*exp(-2*pi*(Umesh.*modkmesh).^2).*cos(2*pi*xyz*HKLt(:,:,S));
            Fhkl(:,S) = sum(F_xyz_hkl,1);
        end
    end
    
    function Fhkl = recalcFhkl(xyz,U)
        Fhkl = zeros(nConv,nSpots);
        for S = 1:nSpots
            [modkmesh, Umesh] = meshgrid(modk(:,S),U);
            B = Batomic_recalc(:,:,S);
            B(freeUiso_recalc,:) = exp(-2*pi*(Umesh.*modkmesh).^2);
            P = Patomic_recalc(:,:,S);
            P(freeAtoms_recalc,:) = cos(2*pi*xyz*HKLt(:,:,S));
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
        
    function Ratios = parList2intensityRatios(PAR)
        Ratios = zeros(1,nSpots);
        XYZfree = reshape(PAR(parIndicesXYZ),[3,nAt_free])';
        Ufree = PAR(parIndicesU);
        Fhkl_free = recalcFhkl(XYZfree,Ufree);
        Fhkl_FED = xExcited*(Fhkl_static + Fhkl_free) + (1-xExcited)*Fhkl_GS;
        Ratios = Fhkl2Icalc(Fhkl_FED)./Igs;
    end

    function X = parList2positions(PAR)
        X = reshape(PAR(parIndicesXYZ),[3,nAt_free])';
    end
end