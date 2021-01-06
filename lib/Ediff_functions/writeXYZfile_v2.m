function output = writeXYZfile_v2(atomSymbols,coords,fileName)
%Line 7 changed from original
nAt = length(atomSymbols);
fileID = fopen(fileName,'a');
fprintf(fileID,'%d \n \n',nAt);
for line = 1:nAt
    atLabel = atomSymbols{line};
    fprintf(fileID,'%s  ',atLabel);
    fprintf(fileID,'%3.5f %3.5f %3.5f \n' ,coords(line,:));
end
fclose(fileID);
end
