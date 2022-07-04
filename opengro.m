function OUT = opengro(filename,nat2read)

fileID = fopen(filename);
if fileID == -1
    disp('.gro file not found')
    return
end
name = string(textscan(fileID,'%s',1));
nAt = str2num(string(textscan(fileID,'%s',1)));
if nargin > 1
    nAt = nat2read;
end
for A = 1:nAt
    grodata{A,:} = textscan(fileID,'%s %s %s %s %s %s %s %s %s', 1);
end
fclose(fileID);
for A = 1:nAt
    labels{A} = string(grodata{A}{2});
    atomTypes{A} = regexp(labels{A},'\D','match');
    labelNumbers{A} = regexp(labels{A},'\d','match');
    atomNumbers(A) = str2num(string(grodata{A}{3})); %should == A
    offset = 0;
    if atomNumbers(A) ~= A
        disp(['Possibile error in atom number ',num2str(atomNumbers(A)),' ~= ',num2str(A)]);
        if length(labelNumbers(A)) >= 6
            offset = -1;
        end
    end
    X(A) = str2num(string(grodata{A}{4+offset}));
    Y(A) = str2num(string(grodata{A}{5+offset}));
    Z(A) = str2num(string(grodata{A}{6+offset}));
    atomNumbers(A) = A;
end

OUT.name = name;
OUT.nAt = nAt;
OUT.labels = labels;
OUT.atomTypes = atomTypes;
OUT.labelNumbers = labelNumbers;
OUT.atomNumbers = atomNumbers;
OUT.XYZ = 10*[X',Y',Z']; %convert from nm to Angstrom
end