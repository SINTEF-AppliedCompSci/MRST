function [x,w,nPts] = fitMoments(x, basis, moments, num, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('reduce', true, 'tol', 1e-12);
    opt = merge_options(opt, varargin{:});

    psi  = basis.psi;
    nDof = basis.nDof;
    nPts = size(x,1);
    
    k       = nPts;
    reduced = true;
%     w = [];
    
%     opts = optimoptions('lsqlin');
%     opts.ConstraintTolerance = 1e-8;
%     opts.OptimalityTolerance = 1e-3;
%     opts.MaxIterations       = 100;
    
     % Matrix of basis functions evalauted at current quadrature points
    p = cell2mat(cellfun(@(p) p(x), psi, 'unif', false));
    [ii, jj] = blockDiagIndex(nPts*ones(num,1), nDof*ones(num,1));
    P = sparse(jj, ii, repmat(p, num,1));

%     I = speye(nPts*num); o = zeros(nPts*num,1); r = Inf(nPts*num,1);
    [w, ~, ~, flag] = lsqnonneg(P, moments);
    flag = flag && all(abs(residual)<opt.tol);
%     [w, ~, ~, flag] = lsqlin(I, o, I, r, P  , moments, [], [], [], opts);
%     [w, ~, ~, flag] = lsqlin(I, o, I, r, P  , moments, o, r, [], opts);
    
    if flag < 0
        error('Moment fitting did not find a solution!');
    end
    
    if opt.reduce
        while k > 0 && reduced 

            % Matrix of basis functions evalauted at current quadrature points
            p = cell2mat(cellfun(@(p) p(x), psi, 'unif', false));
            s = sum(reshape(p.^2, nPts, nDof), 2);
            [~, ix] = sort(s);
            reduced = false;
            xPrev = x;
            wPrev = w;

            for pNo = 1:numel(ix)
                x(ix(pNo),:) = [];
                nPts = size(x,1);
                p = cell2mat(cellfun(@(p) p(x), psi, 'unif', false));
                [ii, jj] = blockDiagIndex(nPts*ones(num,1), nDof*ones(num,1));
                P = sparse(jj, ii, repmat(p, num,1));

%                 I = speye(nPts*num); o = zeros(nPts*num,1); r = Inf(nPts*num,1);
%                 [w, ~, ~, flag] = lsqlin(I, o, I, r, P  , moments, [], [], [], opts);
                [w, ~, ~, flag] = lsqnonneg(P, moments);
                flag = flag && all(abs(residual)<opt.tol);
%                 [w, ~, ~, flag] = lsqlin(I, o, I, r, P  , moments, o, r, [], opts);

                if flag > 0
                    k       = k-1;
                    reduced = true;
                    break
                else
                    x = xPrev;
                    w = wPrev;
                end    
            end

        end
    end
    
    nPts = size(x,1);
    
end
