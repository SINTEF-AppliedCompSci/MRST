function state = initCompositionalState(model, p, T, s0, z0, eos)
% Initialize a compositional state given initial composition
%
% SYNOPSIS:
%   state = initCompositionalState(model, p, T, s0, z0)
%
% PARAMETERS:
%   model - Compositional model, or the reservoir grid itself.
%   p     - Pressures as a column vector
%   T     - Temperatures as a column vector
%   s0    - Initial saturation. Any compositional phases will be flashed, so
%           in that case only the sum matters. Only required for more than
%           two phases present.
%   z     - Composition as a matrix with number of rows equal to the number
%           of components.
%   eos   - EquationOfStateModel instance if model is not the first input
%           argument (legacy syntax). Otherwise optional.
%
% RETURNS:
%   state - Intialized state (after flash)
%

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
    hasModel = ~isstruct(model);
    if hasModel
        % First argument is model
        assert(isa(model, 'ThreePhaseCompositionalModel'));
        if nargin == 5
            eos = model.EOSModel;
        end
        G = model.G;
    else
        G = model;
        assert(nargin > 5);
    end
    state = initResSol(G, p, 1);
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
    sz = size(s0, 2);
    if size(s0, 1) == 1
        s0 = repmat(s0, G.cells.num, 1);
    end
    if hasModel
        % We know the model and can then figure out what phase goes where.
        nph = model.getNumberOfPhases();
        ix = model.getEoSPhaseIndices();
        if nph > 2
            if sz == nph - 2
                % We have value for the non-EoS-phases
                fill = s0;
                s0 = zeros(G.cells.num, nph);
            elseif sz == nph
                % We have one value per phase
                fill = sum(s0, 2) - sum(s0(:, ix), 2);
            else
                error('Bad!');
            end
        else
            fill = 0;
            s0 = zeros(G.cells.num, nph);
        end
    elseif sz == 0
        % Assume two-phase since we got no saturation and no model
        fill = 0;
        s0 = zeros(G.cells.num, 2);
    else
        % We do not know what kind of model we have, just assume that the
        % water phase is the non-EoS
        assert(sz > 1, 'Must have multiple columns in saturations');
        ix = [1, 2] + double(sz > 2);
        fill = sum(s0, 2) - sum(s0(:, ix), 2);
    end
    s0(:, ix) = bsxfun(@times, 1 - fill, [sL, sV]);
    state.s = s0;
end