function state = initCompositionalState(G, p, T, s0, z0, eos)
% Initialize a compositional state given initial composition
%
% SYNOPSIS:
%   state = initCompositionalState(G, p, T, s0, z0, eos)
%
% PARAMETERS:
%   G   - Grid for which the state is to be constructed.
%   p   - Pressures as a column vector
%   T   - Temperatures as a column vector
%   s0  - Initial saturation. Any compositional phases will be flashed, so
%         in that case only the sum matters.
%   z   - Composition as a matrix with number of rows equal to the number
%         of components.
%   eos - EquationOfState derived class instance.
%
% RETURNS:
%   state - Intialized state (after flash)
%

%{
Copyright 2009-2017 SINTEF Digital, Mathematics & Cybernetics.

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

    state = initResSol(G, p, s0);
    state.T = repmat(T, G.cells.num, 1);
    if size(z0, 1) == G.cells.num
        state.components = z0;
    else
        state.components = repmat(z0, G.cells.num, 1);
    end
    nls = getDefaultFlashNonLinearSolver();
    state = eos.validateState(state);
    [state, report] = nls.solveTimestep(state, 1000*year, eos);
    if ~report.StepReports{1}.Converged
        state = eos.updateAfterConvergence(state0, state, dt, struct());
    end
    [sL, sV] = eos.computeSaturations(nan, nan, state.x, state.y, state.L, state.Z_L, state.Z_V);
    sHydrocarbon = [sL, sV];
    if size(state.s, 2) == 3
        % Water is present
        sW = state.s(:, 1);
        state.s = [sW, bsxfun(@times, sHydrocarbon, 1-sW)];
    else
        state.s = sHydrocarbon;
    end
    
end