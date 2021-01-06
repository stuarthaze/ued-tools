function Iratios = fCalcR_Uiso(P,sizeX,sizeY,HKL,modKhkl,Sym,atomicNumbers,F_at,Wts,xExcited,indicies2refine,Fgs_Uiso,Igs,Fes_fixed_Uiso)
addpath('\\win.desy.de\home\hayess\My Documents\MATLAB\lib\Model_Refinement');
[X, U] = fSeparateParameters(P,sizeX,sizeY);
Iratios = fCalcI_ES_Uiso(X,U)./Igs;
    function Intensities = fCalcI_ES_Uiso(x,u)
        atomicNumbersFree = atomicNumbers(indicies2refine);
        F_ES = xExcited*fCalcFhkl_invSym_Uiso(atomicNumbersFree,x(indicies2refine,:),HKL,F_at(indicies2refine,:,:),Sym,u(indicies2refine),modKhkl)+...
            xExcited*Fes_fixed_Uiso + (1-xExcited)*Fgs_Uiso;
        Intensities = sum(abs(F_ES).^2.*Wts,1);
    end
end
 