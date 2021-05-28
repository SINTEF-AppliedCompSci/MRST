classdef OverallCompositionCompositionalModel < ThreePhaseCompositionalModel
    % Overall composition model for compositional problems
    %
    % SYNOPSIS:
    %   model = OverallCompositionCompositionalModel(G, rock, fluid, compFluid)
    %
    % DESCRIPTION:
    %   The overall composition model relies on primary variables
    %   pressure and overall compositions for all components. As a
    %   consequence, the the phase behavior and mixing is offloaded to the
    %   equation of state directly. A requirement for this model is that
    %   the equation of state is fullfilled at each Newton iteration.
    %
    % PARAMETERS:
    %   G         - Grid structure
    %   rock      - Rock structure for the reservoir
    %   fluid     - The flow fluid, containing relative permeabilities,
    %               surface densities and flow properties for the
    %               aqueous/water phase (if present)
    %   compFluid - CompositionalMixture instance describing the species
    %               present.
    %
    % RETURNS:
    %  model - Initialized class instance
    %
    % SEE ALSO:
    %   `ThreePhaseCompositionalModel`, `NaturalVariablesCompositionalModel`

    properties
        
    end
    
    methods
        function model = OverallCompositionCompositionalModel(varargin)
            model = model@ThreePhaseCompositionalModel(varargin{:});
        end
        
        
        function [xM, yM, sL, sV, rhoL, rhoV, muL, muV, report] = computeTwoPhaseFlowProps(model, state, p, temp, z)
            isADI = isa(p, 'ADI') || isa(temp, 'ADI') || any(cellfun(@(x) isa(x, 'ADI'), z));
            report = struct();
            if isADI
                t1 = tic();
                [x, y, L] = model.EOSModel.getPhaseFractionAsADI(state, p, temp, z);
                report.t_derivatives = toc(t1);
                t2 = tic();
                [Z_L, Z_V] = model.EOSModel.getCompressibilityAndFugacity(p, temp, x, y, z, [], []);
                report.t_compressibility = toc(t2);
            else
                [x, y, L, Z_L, Z_V] = model.getProps(state, 'x', 'y', 'L', 'Z_L', 'Z_V');
                x = ensureMinimumFraction(x);
                y = ensureMinimumFraction(y);
                x = expandMatrixToCell(x);
                y = expandMatrixToCell(y);
                [report.t_derivatives, report.t_compressibility] = deal(0);
            end
            
            t1 = tic();
            eos = model.EOSModel;
            xM = eos.getMassFraction(x);
            yM = eos.getMassFraction(y);
            report.t_massfraction = toc(t1);
            
            t2 = tic();
            rhoL = model.PropertyModel.computeDensity(eos, p, x, Z_L, temp, true);
            rhoV = model.PropertyModel.computeDensity(eos, p, y, Z_V, temp, false);
            report.t_density = toc(t2);
            
            [sL, sV] = model.EOSModel.computeSaturations(rhoL, rhoV, x, y, L, Z_L, Z_V);
            
            t3 = tic();
            if nargout > 6
                muL = model.PropertyModel.computeViscosity(eos, p, x, Z_L, temp, true);
                muV = model.PropertyModel.computeViscosity(eos, p, y, Z_V, temp, false);
            end
            report.t_viscosity = toc(t3);
        end
        
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsCompositional(state0, state, model, dt, ...
                            drivingForces, varargin{:});
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            state0 = state;
            
            var0 = problem.primaryVariables;
            vars = var0;
            removed = false(size(vars));
            s_hc = 1;
            extraPhases = model.getNonEoSPhaseNames();
            for i = 1:numel(extraPhases)
                sn = ['s', extraPhases(i)];
                phaseIx = strcmpi(vars, sn);
                if any(phaseIx)
                    state = model.updateStateFromIncrement(state, dx{phaseIx}, problem, sn, inf, model.dsMaxAbs);
                    state = model.capProperty(state, sn, 0, 1);
                    removed(phaseIx) = true;
                    vars = vars(~phaseIx);
                end
                s_hc = s_hc - model.getProp(state, sn);
            end
            % Components
            cnames = model.EOSModel.getComponentNames();
            ncomp = numel(cnames);
            ok = false(ncomp, 1);

            z = state.components;
            rm = 0;

            for i = 1:ncomp
                name = lower(cnames{i});
                cix = strcmpi(var0, name);
                if any(cix)
                    z0 = z(:, i);

                    dz = dx{cix};
                    if isfinite(model.dzMaxAbs)
                        dz = sign(dz).*min(abs(dz), model.dzMaxAbs);
                    end
                    z(:, i) = min(max(z0 + dz, 0), 1);

                    ok(i) = true;
                    [vars, ix] = model.stripVars(vars, {name});
                    removed(~removed) = removed(~removed) | ix;

                    rm = rm - (z(:, i) - z0);
                end
            end
            if any(ok)
                % We had components as active variables somehow
                assert(nnz(~ok) == 1)
                z(:, ~ok) = min(max(z(:, ~ok) + rm, 0), 1);
                z = bsxfun(@rdivide, z, sum(z, 2));
                state.components = z;
                dz = model.computeChange(state.components - state0.components, s_hc);
                state.dz = dz;
            else
                state.dz = zeros(1, ncomp);
            end

            % Parent class handles almost everything for us
            problem.primaryVariables = vars;
            dx(removed) = [];
            [state, report] = updateState@ThreePhaseCompositionalModel(model, state, problem, dx, drivingForces);
            
            if problem.iterationNo == 1
                state.switched = false(model.G.cells.num, 1);
                state.switchCount = zeros(model.G.cells.num, 1);
            end
            twoPhase0 = state.L < 1 & state.L > 0;
            % Set minimum overall composition
            minz = model.EOSModel.minimumComposition;
            state.components = ensureMinimumFraction(state.components, minz);
            % Update saturations etc using flash routine
            state = model.computeFlash(state, problem.dt, problem.iterationNo);
            % Set increments in phase compositions as well
            dx = model.computeChange(state0.x - state.x, s_hc);
            dy = model.computeChange(state0.y - state.y, s_hc);
            dxy = max(dx, dy);
            state.dz = max(state.dz, dxy);
            
            twoPhase = state.L < 1 & state.L > 0;
            switched = twoPhase0 ~= twoPhase;
            
            dispif(model.verbose > 1, '%d pure liquid, %d pure vapor, %d two-phase\n', nnz(state.L == 0), nnz(state.L == 1), nnz(twoPhase));
            state.switchCount = state.switchCount + double(switched);
            % Set minimum phase composition
            state.x = ensureMinimumFraction(state.x, minz);
            state.y = ensureMinimumFraction(state.y, minz);
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
