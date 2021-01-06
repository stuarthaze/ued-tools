function output = writeXYZfile(atomSymbols,coords,fileName)
nAt = length(atomSymbols);
fileID = fopen(fileName,'a');
fprintf(fileID,'%d\n\n',nAt);
for line = 1:nAt
    atLabel = atomSymbols{line}{:};
    fprintf(fileID,'%s  ',atLabel);
    fprintf(fileID,'%3.5f %3.5f %3.5f\n' ,coords(line,:));
end
fclose(fileID);
end
