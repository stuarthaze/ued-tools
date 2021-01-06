function DistancesMat = saveCoordinates2view(CRYSTAL, Refinement, k0, Dists2check)

timeDelays = Refinement.timeDelays;
tIndices = 1:length(timeDelays);
XYZ_T = Refinement.XYZ_T;
axes = CRYSTAL.axes;
atomTypes = CRYSTAL.atomTypes;
nAt = CRYSTAL.nAt;
nT = length(timeDelays);

fCalc_dR = @(X) sqrt(sum(X.^2,2))';
% indiciesBent = 1:2;
% indiciesStr = 3:4;
nDists = size(Dists2check,1);
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
    for at1 = 1:nAt
        x1 = XYZ_cartesian(at1,:);
        for at2 = 1:nAt
            x2 = XYZ_cartesian(at2,:);
            dists(at1,at2,T) = fCalc_dR(x1-x2);
        end
    end
    for d = 1:nDists
        DistancesMat(T,d) = dists(Dists2check(d,1),Dists2check(d,2),T);
    end 
end
% % Calculate values for ground state (GS)
% for at1 = 1:nAt
%     for at2 = 1:nAt
%         dists_GS(at1,at2) = fCalc_dR(cartCoords_gs(at1,:)-cartCoords_gs(at2,:));
%     end
% end
%

figure();
hold on;
plot(timeDelays(tIndices), DistancesMat(:,1),'.-b');
plot(timeDelays(tIndices), DistancesMat(:,2),'x-b');

plot(timeDelays(tIndices), DistancesMat(:,3),'.-k');
plot(timeDelays(tIndices), DistancesMat(:,4),'x-k');

title('Iodine-Iodine distances');
legend('bent 1-2', 'bent 1-3', 'str 4-5','str 4-6');

%
save('dists.txt', 'DistancesMat', '-ascii');
%
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