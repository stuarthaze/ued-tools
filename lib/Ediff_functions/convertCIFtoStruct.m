function CrystalStruct = convertCIFtoStruct(fileDirName)

addpath('\\win.desy.de\home\hayess\My Documents\MATLAB\lib\Crystallography_fileIO');
cifData = cif(fileDirName);
XYZ(:,1) = cifData.atom_site_fract_x;
XYZ(:,2) = cifData.atom_site_fract_y;
XYZ(:,3) = cifData.atom_site_fract_z;
CrystalStruct.XYZ = XYZ;

atomTypes = cifData.atom_site_type_symbol;
CrystalStruct.atomTypes = atomTypes;
CrystalStruct.atomicNumbers = atomicSymbols2Numbers(atomTypes);
%----------------
deg2rad = pi/180.0;
CrystalStruct.a = cifData.cell_length_a;
CrystalStruct.b = cifData.cell_length_b;
CrystalStruct.c = cifData.cell_length_c;
CrystalStruct.alpha_deg = cifData.cell_angle_alpha;
CrystalStruct.beta_deg = cifData.cell_angle_beta;
CrystalStruct.gamma_deg = cifData.cell_angle_gamma;
CrystalStruct.V = cifData.cell_volume;

CrystalStruct.alpha = deg2rad*cifData.cell_angle_alpha;
CrystalStruct.beta = deg2rad*cifData.cell_angle_beta;
CrystalStruct.gamma = deg2rad*cifData.cell_angle_gamma;

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
CrystalStruct.XYZ_cif_cart = fractionals2cartesians(XYZ,axes);
CrystalStruct.nAt = size(CrystalStruct.XYZ,1);
%------------------

% Crystal symmetry - 
% For TBAT, only P-1 symmetry, so refinement needs only identity operation 
CrystalStruct.Sym.point = [1, 1, 1];
CrystalStruct.Sym.trans = [0, 0, 0];
CrystalStruct.nSym = size(CrystalStruct.Sym.point,1);
% To generate full structure, use also inversion:
SymFull.point = [[1, 1, 1];[-1, -1, -1]];
SymFull.trans =[[0,0,0];[0,0,0]];
nSymFull = size(SymFull.point,1);


CrystalStruct.Uiso = cifData.atom_site_U_iso_or_equiv;

%----------------------------------
% Anisotropic displacement params
%----------------------------------
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