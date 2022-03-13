function basis = dgBasis(dim, degree, type, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    opt = struct('k', []);
    opt = merge_options(opt, varargin{:});
    
    sz = size(degree);
    assert(all(sz == 1) || sz(2) == dim);
    
    perDim    = numel(degree) > 1;

    if isempty(opt.k)
        maxDegree = max(degree);
        k = zeros(polyDim(maxDegree, 1).^dim, dim);
        for dNo = 1:dim   
            kk = repmat(0:maxDegree, polyDim(maxDegree, 1).^(dNo-1), ...
                      polyDim(maxDegree, 1).^(dim - dNo));
            k(:,dNo) = kk(:);
        end
        k = k(sum(k,2) <= maxDegree,:);
    else
        k = opt.k;
        maxDegree = max(sum(k,2));
    end
    [~, ix] = sort(sum(k,2));
    k       = k(ix,:);
    
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
    
    basis = struct('psi'    , {psi}    , ...
                   'gradPsi', {gradPsi}, ...
                   'degree' , maxDegree, ...
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
