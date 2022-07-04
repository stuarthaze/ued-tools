function Imol = calcImol_v3(srange,Fs,coord,coherLength)
% S.Hayes 8.June.2020
% Factor of 2 included
% Added: Coherence length option
ns = length(srange);
nat = length(coord(:,1));
Imol = zeros(1,ns);
testdim = size(srange);
if testdim(1) > testdim(2)
    srange = srange';
end

if nargin == 4 && ~isempty(coherLength)
    useCoherence = 1;
    coherFact = -1/(coherLength^2);
else
    useCoherence = 0;
end

for a = 1:(nat-1)
    for b = (a+1):nat
        dx = coord(a,:) - coord(b,:);
        r2rd = sum(dx.^2);
        rab = sqrt(r2rd);
        if useCoherence 
            Imolab = exp(r2rd*coherFact)*Fs(a,:).*Fs(b,:).*sin(srange*rab)./(srange*rab);
        else
            Imolab = Fs(a,:).*Fs(b,:).*sin(srange*rab)./(srange*rab);
        end
        
        if srange(1) == 0 && useCoherence
            Imolab(1) = exp(r2rd*coherFact)*Fs(a,1).*Fs(b,1);
        elseif srange(1) == 0
            Imolab(1) = Fs(a,1).*Fs(b,1);
        end
        Imol = Imol + 2*Imolab;
    end
end
end
