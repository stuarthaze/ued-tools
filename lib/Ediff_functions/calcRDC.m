function [rdc rrange] = calcRDC(sMs,srange,dr,rmax,damping)
   rrange = 0:dr:rmax;
   mRDC = zeros(length(srange),length(rrange));
   for ir = 1:length(rrange)
       r = rrange(ir);
       mRDC(:,ir) = sMs.*sin(r*srange).*exp(-damping*srange.^2);
   end
   rdc = sum(mRDC);
end
