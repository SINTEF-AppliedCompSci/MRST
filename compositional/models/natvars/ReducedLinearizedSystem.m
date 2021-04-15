classdef ReducedLinearizedSystem < LinearizedProblem
    % Linearized system which performs a Schur complement before solving
    % the linear system
    properties
        keepNum % Number of equations to keep after elimination
        reorder % Reordering indices for the variables before elimination
        B % Kept upper block
        C % Eliminated upper block
        D % First lower block
        E % Second lower block
        f % First part of right-hand-side
        h % Last part of right-hand-side
        E_U % Upper part of LU-factorized E-matrix
        E_L % Lower part of LU-factorized E-matrix
    end
    methods
        function problem = ReducedLinearizedSystem(varargin)
            problem = problem@LinearizedProblem(varargin{:});
        end
        
        function p = assembleSystem(p)
           preassembly = ~isempty(p.A);
           p = assembleSystem@LinearizedProblem(p);
           if ~preassembly
               A = p.A;
               b = p.b;
               if ~isempty(p.reorder)
                   A = A(:, p.reorder);
               end
               start = 1:p.keepNum;
               if isempty(A)
                   return
               end
               if 1
                   [ix, jx, vx] = find(A);
                   if any(~isfinite(vx))
                       warning('Non-finite values in matrix before Schur-complement reduction.');
                   end
                   if any(~isfinite(b))
                       warning('Non-finite values in right-hand side before Schur-complement reduction.');
                   end
                   n = size(A, 2);
                   keep = false(n, 1);
                   keep(start) = true;
                   nk = p.keepNum;

                   keepRow = keep(ix);
                   keepCol = keep(jx);
                   kb = keepRow & keepCol;
                   p.B = sparse(ix(kb), jx(kb), vx(kb), nk, nk);

                   kc = keepRow & ~keepCol;
                   p.C = sparse(ix(kc), jx(kc) - nk, vx(kc), nk, n - nk);

                   kd = ~keepRow & keepCol;
                   p.D = sparse(ix(kd) - nk, jx(kd), vx(kd), n - nk, nk);

                   ke = ~keepRow & ~keepCol;
                   p.E = sparse(ix(ke) - nk, jx(ke) - nk, vx(ke), n - nk, n - nk);
                   p.f = b(keep);
                   p.h = b(~keep);
               else
                   stop = (p.keepNum+1):size(A, 2);
                   p.B = A(start, start);
                   p.C = A(start, stop);
                   p.D = A(stop, start);
                   p.E = A(stop, stop);

                   p.f = b(start);
                   p.h = b(stop);
               end
               [L, U] = lu(p.E);
               p.A = p.B - p.C*(U\(L\p.D));
               p.b = p.f - p.C*(U\(L\p.h));
               p.E_U = U;
               p.E_L = L;
           end
        end
        
        function [p, report] = processResultAfterSolve(problem, p, report)
           if isempty(problem.E)
               return
           end
           % s = problem.E\(problem.h - problem.D*p);
           s = problem.E_U\(problem.E_L\(problem.h - problem.D*p));
           if any(~isfinite(s))
               warning('Recovered values were non-finite.');
           end
           p = [p; s];
           if ~isempty(problem.reorder)
                backmap = zeros(numel(p), 1);
                backmap(problem.reorder) = (1:numel(p))';
                p = p(backmap);
           end
        end
    end
end

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
