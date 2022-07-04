function Imol = calcImol2(srange,Fs,coord)
ns = length(srange);
nat = length(coord(:,1));
Imol = zeros(1,ns);
testdim = size(srange);
if testdim(1) > testdim(2)
    srange = srange';
end
for a = 2:nat
    for b = 1:(a-1)
        rab2 = 0.0;
        for i = 1:3
            rab2 = rab2 + (coord(a,i) - coord(b,i))^2;
        end
        rab = sqrt(rab2);
        Imolab = Fs(a,:).*Fs(b,:).*sin(srange*rab)./(srange*rab);
        Imol = Imol + Imolab;
    end
end
end
