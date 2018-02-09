function [k, n] = dgBasis(degree, dim)
    
    if degree < 0
        n = 0;
        k = [];
    else
        n = nchoosek(degree + dim, degree);
        
        k = zeros(polyDim(degree, 1).^dim, dim);
        for dNo = 1:dim   
            kk = repmat(0:degree, polyDim(degree, 1).^(dNo-1)    , ...
                                  polyDim(degree, 1).^(dim - dNo));
            k(:,dNo) = kk(:);
        end
        k = k(sum(k,2) <= degree,:);
    end
    
end

function n = polyDim(degree, dim)

    if degree < 0
        n = 0;
    else
        n = nchoosek(degree + dim, degree);
    end
    
end
        