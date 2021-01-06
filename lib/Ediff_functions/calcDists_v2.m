function Rab = calcDists_v2(AB_list,XYZ)
    A = AB_list(:,1);
    B = AB_list(:,2);
    Rab = calcR(XYZ(A,:)-XYZ(B,:));
    function R = calcR(XYZ)
        R = sqrt(sum(XYZ.^2,2));
    end
end