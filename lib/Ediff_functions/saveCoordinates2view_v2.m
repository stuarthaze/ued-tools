function OUT = saveCoordinates2view_v2(CRYSTAL, Refinement, k0)

timeDelays = Refinement.timeDelays;
tIndices = 1:length(timeDelays);
XYZ_T = Refinement.XYZ_T;
axes = CRYSTAL.axes;
atomTypes = CRYSTAL.atomTypes;
nAt = CRYSTAL.nAt;
nT = length(timeDelays);

%------------------------------------------
k0_Cartesian = fractionals2cartesians(k0,axes);
V_vert = [1,0,0];
XYZ_GS_Cartesian = fractionals2cartesians(CRYSTAL.XYZ,axes);
[XYZ_GS_cartesian_k0view, Tmat] = transformCoords_fromViewVector(XYZ_GS_Cartesian,k0_Cartesian,V_vert);
writeXYZfile_v2(atomTypes,XYZ_GS_cartesian_k0view,'output_GS_cartesians_k0view.xyz');

for T = tIndices
    T
    XYZ_cartesian = fractionals2cartesians(XYZ_T(1:nAt,:,T),axes);
    writeXYZfile_v2(atomTypes,XYZ_cartesian,'output_cartesians.xyz');
	[XYZ_cartesian_k0_view, Tmat] = transformCoords_fromViewVector(XYZ_cartesian,k0_Cartesian,V_vert);
	writeXYZfile_v2(atomTypes,XYZ_cartesian_k0_view,'output_cartesians_k0.xyz');
end


% Add ground-state unit cells at +- a-axis
SymStruct_311super.point = ones(3);
SymStruct_311super.trans = zeros(3);
SymStruct_311super.trans(2,1) = -1;
SymStruct_311super.trans(3,1) = 1;
%[XYZ_311_supercell, atomNumbers_311_supercell] = generateEquivalentAtoms(XYZ_cif, atomicNumbers, SymStruct_311super);
[XYZ_311_supercell, atomTypes_311_supercell] = generateEquivalentAtoms(CRYSTAL.XYZ, atomTypes, SymStruct_311super);
% Save supercell coords
XYZ_311_supercell_cart = fractionals2cartesians(XYZ_311_supercell, axes);
XYZ_311_supercell_k0view = transformCoords_fromViewVector(XYZ_311_supercell_cart,k0_Cartesian,V_vert);
writeXYZfile_v2(atomTypes_311_supercell, XYZ_311_supercell_k0view, 'XYZ_311_super_k0view.xyz');


XYZ_T_supercell = zeros([nAt*3,3,nT]);

for T = tIndices
    XYZ_T_supercell(:,:,T) = XYZ_311_supercell;
    XYZ_T_supercell(1:nAt,:,T) = XYZ_T(:,:,T);
    XYZ_T_super_cart(:,:,T) = fractionals2cartesians(XYZ_T_supercell(:,:,T), axes);
    XYZ_T_super_k0view(:,:,T) = transformCoords_fromViewVector(XYZ_T_super_cart(:,:,T),k0_Cartesian,V_vert);
    writeXYZfile_v2(atomTypes_311_supercell, XYZ_T_super_k0view(:,:,T), 'XYZ_T_super_k0view.xyz');
end
end