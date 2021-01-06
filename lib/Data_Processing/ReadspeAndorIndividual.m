function [dat] = ReadspeAndorIndividual(filename)
% Individual backup frames are saved as 'uint16' format
fid=fopen(filename); 
Im1=fread(fid,'uint16'); % reads the entire file, including metadata

px=Im1(2);
py=Im1(3);

Im=Im1(4:end);
Z=reshape(Im,px,py);
dat=double(Z);

fclose(fid);
end