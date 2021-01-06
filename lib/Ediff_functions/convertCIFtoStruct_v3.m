function CrystalStruct = convertCIFtoStruct_v3(fileDirName)

addpath('\\win.desy.de\home\hayess\My Documents\MATLAB\lib\Crystallography_fileIO');
cifData = cif(fileDirName);
UVW(:,1) = cifData.atom_site_fract_x;
UVW(:,2) = cifData.atom_site_fract_y;
UVW(:,3) = cifData.atom_site_fract_z;
CrystalStruct.UVW_cif = UVW;
CrystalStruct.UVW = UVW;

atomTypes = cifData.atom_site_type_symbol;
CrystalStruct.atomTypes = atomTypes;
CrystalStruct.atomicNumbers = atomicSymbols2Numbers(atomTypes);
if min(CrystalStruct.atomicNumbers < 1)
    disp('Error: atom not found');
end
%----------------
deg2rad = pi/180.0;
a = cifData.cell_length_a;
CrystalStruct.a = a;
b = cifData.cell_length_b;
CrystalStruct.b = b;
c = cifData.cell_length_c;
CrystalStruct.c = c;
CrystalStruct.alpha_deg = cifData.cell_angle_alpha;
CrystalStruct.beta_deg = cifData.cell_angle_beta;
CrystalStruct.gamma_deg = cifData.cell_angle_gamma;

CrystalStruct.alpha = deg2rad*cifData.cell_angle_alpha;
alpha = CrystalStruct.alpha;
CrystalStruct.beta = deg2rad*cifData.cell_angle_beta;
beta = CrystalStruct.beta;
CrystalStruct.gamma = deg2rad*cifData.cell_angle_gamma;
gamma = CrystalStruct.gamma;

if cifData.cell_volume
    CrystalStruct.V = cifData.cell_volume;
else
    CrystalStruct.V = a*b*c*sqrt(1 + 2*cos(alpha)*cos(beta)*cos(gamma) - cos(alpha)^2 - cos(beta)^2 - cos(gamma)^2);
end

% Calculate axes in Cartesian basis
axes.a = CrystalStruct.a*[1 0 0];
axes.b = CrystalStruct.b*[cos(CrystalStruct.gamma), sin(CrystalStruct.gamma), 0];
C = CrystalStruct.c;
axes.c = [C*cos(CrystalStruct.beta), C*(cos(CrystalStruct.alpha)-cos(CrystalStruct.beta)*cos(CrystalStruct.gamma))/sin(CrystalStruct.gamma), ...
    CrystalStruct.V/(CrystalStruct.a*CrystalStruct.b*sin(CrystalStruct.gamma))];
CrystalStruct.vol = dot(axes.a,cross(axes.b,axes.c));
% %-----------------------------
% % reciprocal space
axes.a_star = cross(axes.b,axes.c)/CrystalStruct.vol;
axes.b_star = cross(axes.c,axes.a)/CrystalStruct.vol;
axes.c_star = cross(axes.a,axes.b)/CrystalStruct.vol;
%------------------
CrystalStruct.axes.a = axes.a;
CrystalStruct.axes.b = axes.b;
CrystalStruct.axes.c = axes.c;
CrystalStruct.axes.a_star = axes.a_star;
CrystalStruct.axes.b_star = axes.b_star;
CrystalStruct.axes.c_star = axes.c_star;

%-------------------
% Transformation matrix (xyz)=(uvw)*T
% NOTE: The order is reversed compared to the usual converntion
% since the xyz coords are a row vector
TR(1,:) = axes.a;
TR(2,:) = axes.b;
TR(3,:) = axes.c;
Tinv = inv(TR);
CrystalStruct.T_FracCart = TR;
CrystalStruct.T_CartFrac = Tinv;

CrystalStruct.XYZ_cif = UVW*TR;
CrystalStruct.XYZ = CrystalStruct.XYZ_cif;
CrystalStruct.nAt = size(CrystalStruct.UVW,1);
%------------------
% For converting hkl -> xyz
T2(1,:) = axes.a_star;
T2(2,:) = axes.b_star;
T2(3,:) = axes.c_star;
CrystalStruct.T_hkl2xyz = T2;
CrystalStruct.T_xyz2hkl = inv(T2);

% Crystal symmetry - 
% For TBAT, only P-1 symmetry, so refinement needs only identity operation 
CrystalStruct.Sym.point = [1, 1, 1];
CrystalStruct.Sym.trans = [0, 0, 0];
CrystalStruct.nSym = size(CrystalStruct.Sym.point,1);
% To generate full structure, use also inversion:
SymFull.point = [[1, 1, 1];[-1, -1, -1]];
SymFull.trans =[[0,0,0];[0,0,0]];
nSymFull = size(SymFull.point,1);

if cifData.atom_site_U_iso_or_equiv
CrystalStruct.Uiso = cifData.atom_site_U_iso_or_equiv;
else
    CrystalStruct.Uiso = zeros([CrystalStruct.nAt,1]);
end

%----------------------------------
% Anisotropic displacement params
%----------------------------------
if cifData.atom_site_aniso_U_11
    Ucif = zeros(CrystalStruct.nAt,6);
    atomIndiciesH = CrystalStruct.atomicNumbers == 1;
    atomIndiciesNonH = not(atomIndiciesH);
    Ucif(atomIndiciesNonH,1) = cifData.atom_site_aniso_U_11;
    Ucif(atomIndiciesNonH,2) = cifData.atom_site_aniso_U_22;
    Ucif(atomIndiciesNonH,3) = cifData.atom_site_aniso_U_33;
    Ucif(atomIndiciesNonH,4) = cifData.atom_site_aniso_U_23;
    Ucif(atomIndiciesNonH,5) = cifData.atom_site_aniso_U_13;
    Ucif(atomIndiciesNonH,6) = cifData.atom_site_aniso_U_12;
    Uijn = fCreateUij_matrix(Ucif);
    CrystalStruct.Uijn = Uijn;
    Uequiv = fUij2Uequiv(Uijn,axes);
end

end