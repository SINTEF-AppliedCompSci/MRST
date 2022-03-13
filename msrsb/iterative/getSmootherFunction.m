function fn = getSmootherFunction(varargin)
%Get function for setting up smoother functions
%
% SYNOPSIS:
%   fn = getSmootherFunction();
%   smoother = fn(A, b);
%
% DESCRIPTION:
%   Generates a function handle on the format @(A, b) which will take in a
%   given matrix and a right hand side and return another function handle
%   on the format @(d, x) which applies a smoother to the system Ax = d.
%
% OPTIONAL PARAMETERS:
%   Type       - Type of smoother. Can be 'jacobi', 'sor', for successive
%                over-relaxation (SOR), or  'ilu' for ilu(0).
%
%   Iterations - How many iterations the smoother should apply. Default 5
%                for Jacobi and SOR, and 1 for ILU(0).
%
%   Omega      - If type = 'sor', parapeter for successive over-relaxation.
%                0 < omega < 2, omega = 1 gives the Gauss-Seidel method.
%
% RETURNS:
%  fn  - Function handle for setting up smoother (see description).
%
%
% SEE ALSO:
%   `incompTPFA`, `solveMultiscaleIteratively`, `twoStepMultiscalePreconditioner`

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


    opt = struct('type', 'jacobi', 'iterations', nan, 'omega', 1);
    opt = merge_options(opt, varargin{:});
    
    switch lower(opt.type)
        case 'jacobi'
            if isnan(opt.iterations)
                opt.iterations = 5;
            end
            fn = @(A, b) getJacobiSmoother(A, b, opt.iterations);
        case {'ilu', 'ilu0'}
            if isnan(opt.iterations)
                opt.iterations = 1;
            end
            fn = @(A, b) getILU0(A, b, opt.iterations);
        case 'sor'
            if isnan(opt.iterations)
                opt.iterations = 5;
            end
            fn = @(A, b) getSOR(A, b, opt.omega, opt.iterations);
        otherwise
            error('Unknown smoother')
    end
end

function smoother = getJacobiSmoother(A, ~, its)
    d = diag(A);
    Ar = A - diag(d);
    % D_inv = diag(1./d);
    n = size(A, 1);
    D_inv = spdiags(1./d, 0, n, n);
    % Apply a single pass of jacobi
    jac = @(d, x) D_inv*(d - Ar*x);
    smoother = @(d) loopfun(d, jac, its);
end

function smoother = getILU0(A, ~, its)
    setup = struct('type', 'nofill');
    
    [L, U] = ilu(A, setup);
    if its == 1
        smoother = @(d) U\(L\d);
    else
        jac = @(d, x) x + U\(L\(d - A*x));
        smoother = @(d) loopfun(d, jac, its);
    end
end

function smoother = getSOR(A, ~, omega, its)
    
    n = size(A, 1);
    U = -triu(A, 1);
    L = -tril(A,-1);
    D = spdiags(diag(A), 0, n, n);
    
    sor = @(d,x) (D - omega*L)\(omega*d + (omega*U + (1-omega)*D)*x);
    smoother = @(d) loopfun(d, sor, its);
    
end

function y = loopfun(d, smooth, its)
    y = zeros(size(d));
    for i = 1:its
        y = smooth(d, y);
    end
end
