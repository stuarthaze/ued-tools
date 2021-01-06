function [xyz2, atoms2, occupancies2] = selectCoordSubset_withOccup(xyz,atoms,occupancies,xLimits,yLimits,zLimits)
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
nAt = length(atoms);
at2indx = 0;
for atm = 1:nAt
    if (x(atm) > xLimits(1)) && (x(atm) < xLimits(2))
        if (y(atm) > yLimits(1)) && (y(atm) < yLimits(2))
            if (z(atm) > zLimits(1)) && (z(atm) < zLimits(2))
                at2indx = at2indx+1;
                xyz2(at2indx,:) = xyz(atm,:);
                atoms2(at2indx) = atoms(atm);
                occupancies2(at2indx) = occupancies(atm);
            end
        end
    end
end
