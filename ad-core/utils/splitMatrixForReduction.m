function sys = splitMatrixForReduction(A, b, n, strategy, doFactor, fullFactor)
% Split matrix A and right-hand side into blocks
%  A = [B, C] b = [f] 
%      [D, E]     [h]

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
    if nargin < 6
        fullFactor = false;
    end
    if nargin < 5
        doFactor = true;
    end
    if isempty(n) || n >= size(b, 1) || isempty(A)
        sys = initStruct();
        return
    end
    if any(~isfinite(b))
        warning('Non-finite values in right-hand side before Schur-complement reduction.');
    end
    
    switch lower(strategy)
        case 'sparse'
            % Perform a sparse reduction
            sys = reduceSparse(A, b, n);
        case 'subset'
            % Use Matlab's subset
            sys = reduceSubset(A, b, n);
        otherwise
            error('Unknown reduction strategy %s', strategy);
    end

    if doFactor
        % Perform LU-factorization of E to make it easier to form the
        % reduced Schur complement
        if fullFactor
            % Call full outputs to emulate backslash (and avoid warnings on
            % Octave.)
            [sys.E_L, sys.E_U, sys.P, sys.Q, sys.D] = lu(sys.E);
        else
            [sys.E_L, sys.E_U] = lu(sys.E);
        end
    end
end

function sys = initStruct()
    sys = struct('B', [], 'C', [], 'D',   [], 'E',   [],...
                 'f', [], 'h', [], 'E_L', [], 'E_U', []);
end


function [ix, jx, vx] = breakdownMatrix(A)
    [ix, jx, vx] = find(A);
    if any(~isfinite(vx))
        warning('Non-finite values in matrix before Schur-complement reduction.');
    end
end

function sys = reduceSparse(A, b, nk)
    sys = initStruct();
    start = 1:nk;
    [ix, jx, vx] = breakdownMatrix(A);
    n = size(A, 2);
    keep = false(n, 1);
    keep(start) = true;
    keepRow = keep(ix);
    keepCol = keep(jx);
    kb = keepRow & keepCol;
    sys.B = sparse(ix(kb), jx(kb), vx(kb), nk, nk);

    kc = keepRow & ~keepCol;
    sys.C = sparse(ix(kc), jx(kc) - nk, vx(kc), nk, n - nk);

    kd = ~keepRow & keepCol;
    sys.D = sparse(ix(kd) - nk, jx(kd), vx(kd), n - nk, nk);

    ke = ~keepRow & ~keepCol;
    sys.E = sparse(ix(ke) - nk, jx(ke) - nk, vx(ke), n - nk, n - nk);
    sys.f = b(keep, :);
    sys.h = b(~keep, :);
end

function sys = reduceSubset(A, b, nk)
    sys = initStruct();
    start = 1:nk;
    stop = (nk+1):size(A, 2);
    sys.B = A(start, start);
    sys.C = A(start, stop);
    sys.D = A(stop, start);
    sys.E = A(stop, stop);

    sys.f = b(start, :);
    sys.h = b(stop, :);
end
