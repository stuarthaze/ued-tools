function R = convolveCircle(A,diameter)

arrayDiam = ceil(diameter);
cen = arrayDiam/2 + 0.5;
Rmat = fGenerate2dRmatrix(cen,cen,[arrayDiam,arrayDiam]);
radius = diameter*0.5;
B = Rmat <= radius;
R = conv2(A,double(B),'same');
end
