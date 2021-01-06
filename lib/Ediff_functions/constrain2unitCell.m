function uvw2 = constrain2unitCell(uvw)
uvw2 = uvw;
negCoords = uvw < 0;
uvw2(negCoords) = uvw(negCoords) + 1;
coordsGT1 = uvw >= 1;
uvw2(coordsGT1) = uvw(coordsGT1) - 1;
end

