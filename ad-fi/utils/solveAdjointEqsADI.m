function [x, ii] = solveAdjointEqsADI(eqs, eqs_p, adjVec, objk, system)
%Undocumented Utility Function

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

numVars = cellfun(@numval, eqs)';
cumVars = cumsum(numVars);
ii = [[1;cumVars(1:end-1)+1], cumVars];

if iscell(objk)
   objk = objk{:};
end
objk = cat(objk);

% Above CAT means '.jac' is a single element cell array.  Extract contents.
rhs  = -objk.jac{1}';
if ~isempty(adjVec)
    % If adjVec is not empty, we are not at the last timestep (first in the
    % adjoint recurrence formulation). This means that we subtract
    % the previous jacobian times the previous lagrange multiplier.
    % Previous here means at time t + 1.
    eqs_p = cat(eqs_p{:});

    % CAT means '.jac' is a single element cell array.
    rhs = rhs - eqs_p.jac{1}'*adjVec;
end
tic
if system.nonlinear.cpr
    [x, its, fl] = cprAdjoint(eqs, rhs, system, 'cprType', system.nonlinear.cprType, 'relTol', ...
                   system.nonlinear.cprRelTol, 'cprEllipticSolver', system.nonlinear.cprEllipticSolver);
else
    eqs = cat(eqs{:});

    % CAT means '.jac' is a single element cell array.
    x = eqs.jac{1}'\rhs;
end
tt = toc;
dispif(false, 'Lin. eq: %6.5f seconds, ', tt);
