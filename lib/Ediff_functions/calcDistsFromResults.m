function distsOut = calcDistsFromResults(distsIn,RESULTS,CRYSTAL,checkStart)

distsOut = distsIn;
XYZT = RESULTS.XYZ_T;
UVWT = RESULTS.UVW_T;
nT = size(XYZT,3);

if nargin >= 4 && checkStart
    for T = 1:nT
        XYZT(RESULTS.atomNumbersRefining,:,T) = RESULTS.dXYZstart(:,:,T) + CRYSTAL.XYZ(RESULTS.atomNumbersRefining,:);
        UVWT(:,:,T) = XYZT(:,:,T)*CRYSTAL.T_CartFrac;
    end
end

A = distsIn.atA;
B = distsIn.atB;

R = zeros([distsIn.num,nT]);
Rgs = zeros([distsIn.num,1]);

for P = 1:distsIn.num
    for T = 1:nT
        xyzA = XYZT(A(P),:,T);
        if distsIn.atB_useES
            uvwB = UVWT(B(P),:,T);
        else
            uvwB = CRYSTAL.UVW(B(P),:);
        end
        uvwB = uvwB + distsIn.atB_displUVW(P,:);
        xyzB = uvwB*CRYSTAL.T_FracCart;
        R(P,T) = norm(xyzA-xyzB);
    end
    % Calculate GS values
    xyzA = CRYSTAL.XYZ(A(P),:);
    uvwB = CRYSTAL.UVW(B(P),:);
    uvwB = uvwB + distsIn.atB_displUVW(P,:);
    xyzB = uvwB*CRYSTAL.T_FracCart;
    Rgs(P) = norm(xyzA-xyzB);
end
distsOut.R = R;

if distsIn.plot
    plot(R');
    title('Inter-molecular distances, ES-GS');
    for P = 1:distsIn.num
        legendtext{P} = num2str(distsIn.atA(P));
    end
    legend(legendtext);
end