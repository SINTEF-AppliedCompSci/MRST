function [k, n, psi, grad_psi] = dgBasis(degree, dim)
    
    if degree < 0
        n = 0;
        k = [];
    else
        n = polyDim(degree, dim);
        
        k = zeros(polyDim(degree, 1).^dim, dim);
        for dNo = 1:dim   
            kk = repmat(0:degree, polyDim(degree, 1).^(dNo-1)    , ...
                                  polyDim(degree, 1).^(dim - dNo));
            k(:,dNo) = kk(:);
        end
        k = k(sum(k,2) <= degree,:);
        
    end
    
    [psi, grad_psi] = deal(cell(n,1));
    for bNo = 1:n
        psi{bNo} = @(x) prod(x.^k(bNo,:),2);
        grad_psi{bNo} = @(x) [k(bNo,1).*prod(x.^[max(k(bNo,1)-1,0), k(bNo,2)         ],2) ...
                              k(bNo,2).*prod(x.^[k(bNo,1)         , max(k(bNo,2)-1,0)],2)];
    end
    
end

function n = polyDim(degree, dim)

    if degree < 0
        n = 0;
    else
        n = nchoosek(degree + dim, degree);
    end
    
end
        