function OUT = saveCoordinates2view_v4(CRYSTAL, Refinement, k0_uvw, Displacements)

timeDelays = Refinement.timeDelays;
tIndices = 1:length(timeDelays);
XYZ_T = Refinement.XYZ_T;
axes = CRYSTAL.axes;
atomTypes = CRYSTAL.atomTypes;
nAt = CRYSTAL.nAt;
nT = length(timeDelays);

if nargin > 3
    nGroups = Displacements.nGroups;
    nAtGrouped = Displacements.nAtGrouped;
    nAtTotal = nGroups*nAtGrouped;
    if ~(nAt == (nAtTotal))
        disp('Warning: Error in number of atoms');
    end
    dUVW_groups = zeros([nAt,3]);
    aa = 0;
    for group = 1:Displacements.nGroups
        for atm = 1:nAtGrouped
            aa = aa + 1;
            dUVW(aa,:) = Displacements.UVW(group,:);
        end
    end
    dXYZ_groups = dUVW*CRYSTAL.T_FracCart;
    for T = tIndices
        XYZ_T(:,:,T) = XYZ_T(:,:,T) + dXYZ_groups;
    end
end

    
    
%------------------------------------------
k0_xyz = k0_uvw*CRYSTAL.T_FracCart;
V_vert = [1,0,0];
XYZ_GS_Cartesian = CRYSTAL.XYZ;
if nargin > 3
    XYZ_GS_Cartesian = XYZ_GS_Cartesian + dXYZ_groups;
end
[XYZ_GS_cartesian_k0view, Tmat] = transformCoords_fromViewVector(XYZ_GS_Cartesian,k0_xyz,V_vert);
writeXYZfile_v2(atomTypes,XYZ_GS_cartesian_k0view,'structures_GS_k0.xyz');

for T = tIndices
    T
    XYZ_cartesian = XYZ_T(1:nAt,:,T);
    writeXYZfile_v2(atomTypes,XYZ_cartesian,'structures_fitted_originalOrientation.xyz');
	[XYZ_cartesian_k0_view, Tmat] = transformCoords_fromViewVector(XYZ_cartesian,k0_xyz,V_vert);
	writeXYZfile_v2(atomTypes,XYZ_cartesian_k0_view,'structures_fitted_k0.xyz');
end


% Add ground-state unit cells at +- a-axis
SymStruct_311super.point = ones(3);
SymStruct_311super.trans = zeros(3);
SymStruct_311super.trans(2,1) = -1;
SymStruct_311super.trans(3,1) = 1;
%[XYZ_311_supercell, atomNumbers_311_supercell] = generateEquivalentAtoms(XYZ_cif, atomicNumbers, SymStruct_311super);
[UVW_311_supercell, atomTypes_311_supercell] = generateEquivalentAtoms(CRYSTAL.UVW, atomTypes, SymStruct_311super);
% Save supercell coords
XYZ_311_supercell_cart = fractionals2cartesians(UVW_311_supercell, axes);
XYZ_311_supercell_k0view = transformCoords_fromViewVector(XYZ_311_supercell_cart,k0_xyz,V_vert);
writeXYZfile_v2(atomTypes_311_supercell, XYZ_311_supercell_k0view, 'XYZ_311_super_k0.xyz');

UVW_T_supercell = zeros([nAt*3,3,nT]);

for T = tIndices
    UVW = XYZ_T(:,:,T)*CRYSTAL.T_CartFrac;
    UVW_T_supercell(:,:,T) = UVW_311_supercell;
    UVW_T_supercell(1:nAt,:,T) = UVW;
    XYZ_T_super_cart(:,:,T) = fractionals2cartesians(UVW_T_supercell(:,:,T), axes);
    XYZ_T_super_k0view(:,:,T) = transformCoords_fromViewVector(XYZ_T_super_cart(:,:,T),k0_xyz,V_vert);
    writeXYZfile_v2(atomTypes_311_supercell, XYZ_T_super_k0view(:,:,T), 'XYZ_T_super_k0.xyz');
end
fclose('all');
end