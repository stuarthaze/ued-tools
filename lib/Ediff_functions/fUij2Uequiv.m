function Ueq = fUij2Uequiv(Uijn,axes)
NumAtms = size(Uijn,3);
Ueq = zeros(NumAtms,1);
A = zeros(3);
A(:,1) = axes.a;
A(:,2) = axes.b;
A(:,3) = axes.c;
Asqrd = A'*A;
norm_star_vec = [norm(axes.a_star),norm(axes.b_star),norm(axes.c_star)];
for N = 1:NumAtms
% first try    Ueq(N) = norm_star_vec*A'*Uijn(:,:,N)*A*norm_star_vec'/3;
    Ueq(N) = norm_star_vec*(Asqrd.*Uijn(:,:,N))*norm_star_vec'/3;
end
end 