function [X,U] = fSeparateParameters(P,Shape1,Shape2)
n1 = Shape1(1)*Shape1(2);
n2 = Shape2(1)*Shape2(2);
X = reshape(P(1:n1),Shape1);
U = reshape(P((n1+1):(n2+n1)),Shape2);
end