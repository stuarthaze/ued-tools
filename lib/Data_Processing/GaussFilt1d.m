function M2 = GaussFilt1d(M,sigma,dim)
% Apply a Gaussian filter to M in dimension dim
ndimM = ndims(M);
sizeM = size(M);
M2 = zeros(sizeM);

if ndimM == 1
    M2 = imgaussfilt(M,sigma);
elseif ndimM == 2
    if dim == 1
        for x = 1:sizeM(2)
            M2(:,x) = imgaussfilt(M(:,x),sigma);
        end
    elseif dim == 2
        for y = 1:sizeM(1)
            M2(y,:) = imgaussfilt(M(y,:),sigma);
        end
    end
    
elseif ndimM == 3
    if dim == 1
        for x = 1:sizeM(2)
            for z = 1:sizeM(3)
                M2(:,x,z) = imgaussfilt(M(:,x,z),sigma);
            end
        end
    elseif dim == 2
        for y = 1:sizeM(1)
            for z = 1:sizeM(3)
                M2(y,:,z) = imgaussfilt(M(y,:,z),sigma);
            end
        end
    elseif dim == 3
        for y = 1:sizeM(1)
            for x = 1:sizeM(2)
                M1d = reshape(M(y,x,:),[1,sizeM(3)]);
                M2(y,x,:) = imgaussfilt(M1d,sigma);
            end
        end
    end
end

end
