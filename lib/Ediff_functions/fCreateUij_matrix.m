function U = fCreateUij_matrix(Uaniso)
% Creates (3,3,NumAtoms) matrix from (NumAtoms,6) independent variables
Natm = size(Uaniso,1);
U = zeros(3,3,Natm);
U(1,1,:) = Uaniso(:,1);
U(2,2,:) = Uaniso(:,2);
U(3,3,:) = Uaniso(:,3);
U(2,3,:) = Uaniso(:,4);
U(1,3,:) = Uaniso(:,5);
U(1,2,:) = Uaniso(:,6);
U(2,1,:) = U(1,2,:);
U(3,1,:) = U(1,3,:);
U(3,2,:) = U(2,3,:);
end