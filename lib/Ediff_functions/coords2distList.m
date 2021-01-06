function distList = coords2dists(coords)
%    dists = zeros(1,nDists);
   nAt = size(coords,1);
   nDists = nAt*(nAt-1)/2;
   distList = zeros(nDists,4);
   count = 0;
   rAB2 = zeros(1,3);
   for a = 1:(nAt-1)
       for b = (a+1):nAt
           count = count + 1;
           distList(count,1:2) = [a, b];
           distList(count, 3) = sqrt(sum((coords(a,:)-coords(b,:)).^2));
       end
   end
end
