function refinedStruct = refineStructure(FED, CRYSTAL, MODEL)

XYZgs = CRYSTAL.XYZ;
XYZgs_cart = fractionals2cartesians(XYZgs,CRYSTAL.axes);
Ugs = CRYSTAL.Uiso;
xExcited = MODEL.xExcited;
HKL = FED.HKL;
modk = fHKL2modK(HKL,CRYSTAL.axes);
nSpots = FED.nSpots;
nConv = FED.nConv;
nTimeDelays = length(FED.timeDelays);

%----- Atomic Scattering Factors  --------------------
% Formatted as [atomIndex,contributionIndex,spotIndex]
Fatomic_all = MODEL.Fatomic;
Fatomic_free = Fatomic_all(MODEL.freeAtoms,:,:);
freeAtoms = MODEL.freeAtoms;
fixedAtoms = ~MODEL.freeAtoms;
Fatomic_fixed = Fatomic_all(fixedAtoms,:,:);

%----- Useful variables ---
nAt_free = sum(MODEL.freeAtoms);
nAt_fixed = CRYSTAL.nAt - nAt_free;
nU_free = sum(MODEL.freeUiso);
freeUiso = MODEL.freeUiso;
fixedUiso = ~freeUiso;
parIndicesXYZ = 1:3*nAt_free;
parIndicesU = (1:nU_free) + parIndicesXYZ(end);

%----- GS structure factors --------
Fhkl_GS = calcFhkl(Fatomic_all,XYZgs,HKL,Ugs);
Fhkl_fixed = calcFhkl(Fatomic_fixed,XYZgs(fixedAtoms,:),HKL,Ugs(fixedAtoms));

Igs = Fhkl2Icalc(Fhkl_GS);

%----- Restraint Functions -----------
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
Pstart = [reshape(XYZgs(MODEL.freeAtoms,:),[1,3*nAt_free]),Ugs(MODEL.freeAtoms)'];
XYZfree_T = zeros(nAt_free,3,nTimeDelays);
XYZfixed_T = zeros(nAt_fixed,3,nTimeDelays);

Ufree_T = zeros(nU_free,nTimeDelays);
options = optimset('Display','iter','MaxFunEvals',1e5,'TolX',1e-6,'MaxIter',1000);
for T = 1:nTimeDelays
    fun2minimize = @(P) [(parList2intensityRatios(P)-FED.intensity_ratios(T,:))./FED.meanUncertainties, fRestraintFunction(P) ];
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

    function Fhkl = calcFhkl(Fat,xyz,hkl,U)
        Fhkl = zeros(size(hkl,1),nSpots);       
        for S = 1:nSpots
            [modkmesh, Umesh] = meshgrid(modk(:,S),U);
            F_xyz_hkl = Fat(:,:,S).*exp(-2*pi*(Umesh.*modkmesh).^2).*cos(2*pi*xyz*hkl(:,:,S)');
            Fhkl(:,S) = sum(F_xyz_hkl,1);
        end
    end
    
    function Icalc = Fhkl2Icalc(Fcalc)
        Icalc = sum(FED.Wts.*Fcalc.^2,1);
    end

    function Ratios = calcRatios(XYZfree,Ufree)
        Ratios = zeros(1,nSpots);
        Fhkl_free = calcFhkl(Fatomic_free,XYZfree,HKL,Ufree);
        Fhkl_FED = xExcited*(Fhkl_fixed + Fhkl_free) + (1-xExcited)*Fhkl_GS;
        Ratios = Fhkl2Icalc(Fhkl_FED)./Igs;
    end
        
    function Ratios = parList2intensityRatios(PAR)
        Ratios = zeros(1,nSpots);
        XYZfree = reshape(PAR(parIndicesXYZ),[nAt_free,3]);
        Ufree = PAR(parIndicesU);
        Fhkl_free = calcFhkl(Fatomic_free,XYZfree,HKL,Ufree);
        Fhkl_FED = xExcited*(Fhkl_fixed + Fhkl_free) + (1-xExcited)*Fhkl_GS;
        Ratios = Fhkl2Icalc(Fhkl_FED)./Igs;
    end

    function X = parList2positions(PAR)
        X = reshape(PAR(parIndicesXYZ),[nAt_free,3]);
    end
end