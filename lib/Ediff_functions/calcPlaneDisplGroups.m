function dXvect = calcPlaneDisplGroups(XYZ,AtomsGroups)
nGroups = length(AtomsGroups);
atmInx = 0;
for grp = 1:nGroups
    nAtGrp = length(AtomsGroups{grp});
    dX{grp} = calcDeviationsFromPlane(XYZ,AtomsGroups{grp});
    indices = (atmInx+1):(atmInx+nAtGrp);
    dXvect(indices) = dX{grp};
    atmInx = atmInx+nAtGrp;
end
end
