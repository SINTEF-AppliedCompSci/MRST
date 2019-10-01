function basis = dgBasis(dim, degree, type)

    perDim = true;
    if numel(degree) == 1
        perDim = false;
        degree = repmat(degree,1,dim);
    end
    assert(numel(degree) == dim);
    maxDegree = max(degree);

    if degree < 0
        nDof = 0;
        k    = [];
    else
        k = zeros(polyDim(maxDegree, 1).^dim, dim);
        for dNo = 1:dim   
            kk = repmat(0:maxDegree, polyDim(maxDegree, 1).^(dNo-1), ...
                      polyDim(maxDegree, 1).^(dim - dNo));
            k(:,dNo) = kk(:);
        end
        k = k(sum(k,2) <= maxDegree,:);
        
        kn = [];
        for d = 0:maxDegree
            kk = k(sum(k,2) == d,:);
            kn = [kn; kk];
        end
        k = kn;
        
        if perDim
            for d = 1:dim
                ix      = k(:,d) > degree(d);
                k(ix,:) = [];
            end
        end
        nDof = size(k,1);
        
        switch type
            case 'poly'
                poly = simplePolynomials(maxDegree);
            case 'legendre'
                poly = legendrePolynomials(maxDegree);
            case 'chebyshev'
                poly = chebyshevPolynomials(maxDegree);                
            otherwise
                error('Unknown basis function class');
        end
        
        [psi, gradPsi] = deal(cell(nDof,1));
        for dofNo = 1:nDof
            p = cell(dim, 1);
            for dNo = 1:dim
                p{dNo} = poly{k(dofNo,dNo)+1};
            end
            psi{dofNo} = tensorProduct(p{:});
            gradPsi{dofNo} = grad(psi{dofNo});
        end

    end
    
    basis = struct('psi'    , {psi}    , ...
                   'gradPsi', {gradPsi}, ...
                   'k'      , k        , ...
                   'nDof'   , nDof     , ...
                   'type'   , type     );
    
end

function n = polyDim(degree, dim)

    if degree < 0
        n = 0;
    else
        n = nchoosek(degree + dim, degree);
    end
    
end