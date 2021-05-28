function [L, x, y, Z_L, Z_V, rhoL, rhoV, reports] = standaloneFlash(p, T, z, EOSModel, varargin)
% Utility for flashing without explicitly forming a state
%
% SYNOPSIS:
%   [L, x, y, Z_L, Z_V] = standaloneFlash(p, T, z, EOSModel)
%
% DESCRIPTION:
%   Wrapper function for solving a EOS flash without dealing with a state.
%
% PARAMETERS:
%   p   - Pressures as a column vector
%   T   - Temperatures as a column vector
%   z   - Composition as a matrix with number of rows equal to the number
%         of components.
%
% SEE ALSO:
%   `EquationOfStateModel`

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
    solver = getDefaultFlashNonLinearSolver(varargin{:});
    
    state = struct();
    [p, T, z] = expandArguments(EOSModel, p, T, z);
    z = max(z, EOSModel.minimumComposition);
    z = bsxfun(@rdivide, z, sum(z, 2));
    state.pressure = p;
    state.T = T;
    state.components = z;
    
    state = EOSModel.validateState(state);
    [state, report] = solver.solveTimestep(state, 1, EOSModel);
    
    L = state.L;
    x = state.x;
    y = state.y;
    Z_L = state.Z_L;
    Z_V = state.Z_V;
    if nargout > 5
        rhoL = EOSModel.PropertyModel.computeDensity(EOSModel, p, x, Z_L, T, true);
        rhoV = EOSModel.PropertyModel.computeDensity(EOSModel, p, y, Z_V, T, false);
    end
    if nargout > 7
        reports = report.StepReports{1}.NonlinearReport;
    end
end

function [p, T, z] = expandArguments(eos, p, T, z)
    p = reshape(p, [], 1);
    T = reshape(T, [], 1);
    ncomp = eos.getNumberOfComponents();
    assert(size(z, 2) == ncomp, ...
        'Component input must have %d columns (the number of components in the eos)', ncomp);
    np = numel(p);
    nt = numel(T);
    nz = size(z, 1);
    
    dims = [np, nt, nz];
    md = max(dims);
    assert(all(dims == 1 | dims == md), 'All inputs must have either 1 row, or equal to the maximum value %d', md);
    if np == 1
        p = repmat(p, md, 1);
    end
    if nt == 1
        T = repmat(T, md, 1);
    end
    if nz == 1
        z = repmat(z, md, 1);
    end
end
