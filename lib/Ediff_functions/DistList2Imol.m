function Imol = DistList2Imol(distances,F,S)
   NS = length(F(1,:));
   nDist = length(distances(:,1));
   Imol = zeros(1,length(S));
   for a = 1:NS
       for b = 1:nDist
           s = S(a);
           F1F2 = F(distances(b,1),a)*F(distances(b,2),a);
           r12 = distances(b,3);
           l12 = distances(b,4);
           I12 = F1F2*exp(-0.5*(l12*s)^2)*sin(r12*s)/(r12*s);
           Imol(a) = Imol(a) + I12;
       end
   end
end
