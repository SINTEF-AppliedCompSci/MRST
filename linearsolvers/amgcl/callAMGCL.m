function [x, err, nIter] = callAMGCL(A, b, varargin)
%Invoke AMGCL Linear Solver Software
%
% DESCRIPTION:
%   Solves the system of simulatenaous linear equations
%
%      Ax = b
%
%   in which the coefficient matrix 'A' is sparse.  Single right-hand side
%   vector.
%
% SYNOPSIS:
%    x              = callAMGCL(A, b)
%    x              = callAMGCL(A, b, 'pn1', pv1, ...)
%   [x, err]        = callAMGCL(...)
%   [x, err, nIter] = callAMGCL(...)
%
% PARAMETERS:
%   A - Coefficient matrix.
%
%   b - System right-hand side vector.
%
% KEYWORD ARGUMENTS:
%   isTransposed  - Whether or not the coefficient matrix is transposed on
%                   input.  This simplifies converting MATLAB's CSC matrix
%                   format to AMGCL's CSR matrix input format.  Default
%                   value: `isTransposed = false`.
%
%   tolerance     - Linear solver tolerance.  Corresponds to parameter
%                   'solver.tol' (relative residual reduction) of AMGCL.
%                   Default value: `tolerance = 1.0e-6`.
%
%   maxIterations - Maximum number of linear iterations.  Overrides
%                   'solver.maxiter' if positive.  Default value:
%                   `maxIterations = 0` (use AMGCL solver's default,
%                   typically 100).
%
%  Additional keyword arguments passed on to function `getAMGCLMexStruct`.
%
% RETURNS:
%   x     - Solution vector.
%
%   err   - Norm of residual at end of solution process.
%
%   nIter - Number of linear iterations.
%
% SEE ALSO:
%   `mldivide`, `amgcl_matlab`, `getAMGCLMexStruct`.

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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

    opt = struct('isTransposed',  false, ...
                 'tolerance',     1e-6, ...
                 'maxIterations',  0);

    [opt, cl_opts] = merge_options(opt, varargin{:});

    amg_opt = getAMGCLMexStruct(cl_opts{:});
    
    % Note the transpose...
    if ~opt.isTransposed
        A = A';
    end

    [x, err, nIter] = ...
       amgcl_matlab(A, b, amg_opt, opt.tolerance, opt.maxIterations, 1);
end
