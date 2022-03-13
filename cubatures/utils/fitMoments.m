function [points,weights,nPts] = fitMoments(x, basis, moments, varargin)
    % Compute cubature by moment fitting, possibly with reduction

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

    opt = struct('reduce', true, 'tol', 1e-10, 'chunkSize', 10, 'equal', false);
    opt = merge_options(opt, varargin{:});

    psi  = basis.psi;
    nDof = basis.nDof;
    n0 = size(x,1);

    ne = numel(moments{1});
    if opt.equal
        ne = 1;
        ePos = 1:2;
    else
        nc = ceil(ne/opt.chunkSize);
        ePos   = round(linspace(0, ne, nc))+1;
    end
    nPts   = zeros(ne,1);
    [points, weights] = deal([]);

    x0     = x;

    for i = 1:numel(ePos)-1
        
        elements = ePos(i):ePos(i+1)-1;
        ne_loc = numel(elements);
        
        fprintf('Compressing quadrature for element %d to %d of %d ... ', ...
                                         min(elements), max(elements), ne);
        tic;
        
        m = cellfun(@(m)m(elements), moments, 'UniformOutput', false);
        tol = 100*eps(mean(m{1}));
        m = vertcat(m{:});
        M = nan(numel(m),1);
        for j = 1:nDof
            M((1:nDof:nDof*ne_loc)+j-1) = m((ne_loc*(j-1) + (1:ne_loc)));
        end
        
        reduced = true;
        x = x0;
        n = n0;
        k = n0;
        w = zeros(n,1);
        
        while k > 0 && reduced 
        
            % Matrix of basis functions evalauted at current quadrature points            
            P = computeBasisMatrix(psi, x, ne_loc, n, nDof);
            [w, ~, residual, flag] = lsqnonneg(P, M);
            
            s = sum(full(P(1:nDof, 1:n)).^2,1);
            if opt.reduce                
                [~, ix] = sort(s);
            else
                ix = [];
            end
            reduced = false;
            xPrev = x;
            wPrev = w;
            nPrev = n;

            for pNo = 1:numel(ix)
                x(ix(pNo),:) = [];
                n = size(x,1);
                
                P = computeBasisMatrix(psi, x, ne_loc, n, nDof);
                [w, ~, residual, flag] = lsqnonneg(P, M);
                flag = flag && all(abs(residual)<tol);
                if flag > 0
                    k       = k-1;
                    reduced = true;
                    break
                else
                    x = xPrev;
                    w = wPrev;
                    n = nPrev;
                end    
            end
            
        end
        
        nPts(elements) = n;
        points  = [points;repmat(x,ne_loc,1)];
        weights = [weights;w];
        
        time = toc;
        fprintf('Compressed from %d to %d points in %f second\n', n0, n, time);
              
    end
    
    if opt.equal
        ne = numel(moments{1});
        nPts    = repmat(n, ne, 1);
        points  = repmat(points, ne, 1);
        weights = repmat(weights, ne, 1);
    end
    
end

function P = computeBasisMatrix(psi, x, ne, n, nDof)
    p = cell2mat(cellfun(@(p) p(x), psi, 'unif', false));
    [ii, jj] = blockDiagIndex(n*ones(ne,1), nDof*ones(ne,1));
    P = sparse(jj, ii, repmat(p, ne,1));
end
