classdef DGBasisFunctions
    
    properties
        type
        psi
        grad_psi
        G
        k
        nDof
    end
    
    methods
        
        function basis = DGBasisFunctions(G, degree, varargin)
            
            basis.type = 'legendre';
            basis.G = G;
            basis = merge_options(basis, varargin{:});
            
            dim = G.griddim;
            nDof = polyDim(degree, G.griddim);
            [psi, grad_psi] = deal(cell(nDof,1));
            if degree < 0
                error('Degree must be non-negative');
            else
                k = zeros(polyDim(degree, 1).^dim, dim);
                for dNo = 1:dim   
                    kk = repmat(0:degree, polyDim(degree, 1).^(dNo-1)    , ...
                    polyDim(degree, 1).^(dim - dNo));
                    k(:,dNo) = kk(:);
                end
                k = k(sum(k,2) <= degree,:);
        
                switch basis.type
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
        
            for dofNo = 1:nDof
                psi{dofNo} = @(x, cells) basis.evaluateBasisFunction(psi{dofNo}, x, cells);
                grad_psi{dofNo} = @(x, cells) basis.evaluateBasisFunctionGradient(grad_psi{dofNo}, x, cells);
            end
                
            basis.psi = psi;
            basis.grad_psi = grad_psi;
            basis.k = k;
            basis.nDof = nDof;
    
        end

        function val = evaluateBasisFunction(basis, p, x, cells)
            xhat = (x - basis.G.cells.centroids(cells,:))./(basis.G.cells.diameters(cells)/(2*sqrt(basis.G.griddim)));
            val  = p(xhat);
        end
        
        function val = evaluateBasisFunctionGradient(basis, p, x, cells)
            G    = basis.G;
            xhat = (x - G.cells.centroids(cells,:))./(G.cells.diameters(cells)/(2*sqrt(G.griddim)));
            val  = p(xhat)./(G.cells.diameters(cells)/(2*sqrt(G.griddim)));
        end
        
    
    end
    
end

function n = polyDim(degree, dim)

    if degree < 0
        n = 0;
    else
        n = nchoosek(degree + dim, degree);
    end
    
end