function basis = dgBasis(degree, dim, type)
    
    nDof = polyDim(degree, dim); 
        [psi, grad_psi] = deal(cell(nDof,1));
    if degree < 0
        nDof = 0;
        k    = [];
    else
        k = zeros(polyDim(degree, 1).^dim, dim);
        for dNo = 1:dim   
            kk = repmat(0:degree, polyDim(degree, 1).^(dNo-1)    , ...
                      polyDim(degree, 1).^(dim - dNo));
            k(:,dNo) = kk(:);
        end
        k = k(sum(k,2) <= degree,:);
        
        switch type
            case 'poly'
                for dofNo = 1:nDof
                    psi{dofNo}      = Polynomial(k(dofNo,:), 1);
                    grad_psi{dofNo} = grad(psi{dofNo});
                end

            case 'legendre'
                leg = legendrePolynomials(degree);
                for dofNo = 1:nDof
                    l = cell(dim, 1);
                    for dNo = 1:dim
                        l{dNo} = leg{k(dofNo,dNo)+1};
                    end
                    psi{dofNo} = combine(l{:});
                    grad_psi{dofNo} = grad(psi{dofNo});
                end
                
            otherwise
                error('Unknown basis function class');
        end

    end
    
    psi_lim = {Polynomial([0,0], 1), ...
               Polynomial([1,0], 1), ...
               Polynomial([0,1], 1)};
    
    basis = struct('psi'         , {psi}     , ...
                   'grad_psi'    , {grad_psi}, ...
                   'psi_lim'     , {psi_lim} , ...
                   'k'           , k         , ...
                   'nDof'        , nDof      , ...
                   'type'        , type      );
    
end

function n = polyDim(degree, dim)

    if degree < 0
        n = 0;
    else
        n = nchoosek(degree + dim, degree);
    end
    
end