classdef BlackOilSolventModel < ThreePhaseBlackOilModel
    % Black-oil solvent model.
    properties
        % Model contains a solvent pseudo-component
        solvent
        % HC residual saturations can be set to the smallest observed over the entire
        % simulaiton history
        hystereticResSat
        % Dynamic end-point scaling should be used in a proper solvent model,
        % but is somewhat unstable. Set to false by default.
        dynamicEndPointScaling
    end

    methods
        function model = BlackOilSolventModel(G, rock, fluid, varargin)
            model = model@ThreePhaseBlackOilModel(G, rock, fluid, varargin{:});

            % Use CNV style convergence 
            model.useCNVConvergence      = true;
            model.hystereticResSat       = false;
            model.dynamicEndPointScaling = false;

            % All phases are present
            model.water   = true;
            model.oil     = true;
            model.gas     = true;
            model.solvent = true;

            model = merge_options(model, varargin{:});

        end

        % --------------------------------------------------------------------%
        function model = validateModel(model, varargin)
            if isempty(model.FacilityModel)
                model.FacilityModel = FacilityModel(model);
            end
            if nargin > 1
                model.FacilityModel = model.FacilityModel.setReservoirModel(model);
                W = varargin{1}.W;
                model.FacilityModel = model.FacilityModel.setupWells(W);
            end
            model = validateModel@ThreePhaseBlackOilModel(model, varargin{:});
            return
        end

        % --------------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
            switch(lower(name))
                case {'solvent', 'ss'}
                    index = 4;
                    fn = 's';
                otherwise
                    % Basic phases are known to the base class
                    [fn, index] = getVariableField@ThreePhaseBlackOilModel(model, name, varargin{:});
            end
        end

        % --------------------------------------------------------------------%
        function vars = getSaturationVarNames(model)
            vars = {'sw', 'so', 'sg', 'ss'};
            ph = model.getActivePhases();
            vars = vars(ph);
        end

        % --------------------------------------------------------------------%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsBlackOilSolvent(state0, state, ...
                    model, dt, drivingForces, varargin{:});

        end

        % --------------------------------------------------------------------%
        function [phNames, longNames] = getPhaseNames(model)
            % Get the active phases in canonical ordering
            tmp = 'WOGS';
            active = model.getActivePhases();
            phNames = tmp(active);
            if nargout > 1
                tmp = {'water', 'oil', 'gas', 'solvent'};
                longNames = tmp(active);
            end
        end

        % --------------------------------------------------------------------%
        function phIndices = getPhaseIndices(model)
            % Get the active phases in canonical ordering
            w = model.water;
            o = model.oil;
            g = model.gas;
            s = model.solvent;
            phIndices = [w, w+o, w+o+g, w+o+g+s];
            phIndices(~model.getActivePhases) = -1;
        end

        % --------------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            vars = problem.primaryVariables;
            removed = false(size(vars));
            if model.disgas || model.vapoil
                % The VO model is a bit complicated, handle this part
                % explicitly.
                state0 = state;

                state = model.updateStateFromIncrement(state, dx, problem, 'pressure', model.dpMaxRel, model.dpMaxAbs);
                state = model.capProperty(state, 'pressure', model.minimumPressure, model.maximumPressure);

                [vars, ix] = model.stripVars(vars, 'pressure');
                removed(~removed) = removed(~removed) | ix;

                % Black oil with dissolution
                so = model.getProp(state, 'so');
                sw = model.getProp(state, 'sw');
                sg = model.getProp(state, 'sg');
                ss = model.getProp(state, 'ss');

                % Magic status flag, see inside for doc
                st = model.getCellStatusVO(state0, so, sw + ss, sg);

                dr = model.getIncrement(dx, problem, 'x');
                dsw = model.getIncrement(dx, problem, 'sw');
                % Interpretation of "gas" phase varies from cell to cell, remove
                % everything that isn't sG updates
                dsg = st{3}.*dr - st{2}.*dsw;
                dss = model.getIncrement(dx, problem, 'ss');

                if model.disgas
                    state = model.updateStateFromIncrement(state, st{1}.*dr, problem, ...
                                                           'rs', model.drsMaxRel, model.drsMaxAbs);
                end

                if model.vapoil
                    state = model.updateStateFromIncrement(state, st{2}.*dr, problem, ...
                                                           'rv', model.drsMaxRel, model.drsMaxAbs);
                end

                dso = -(dsg + dsw + dss);
                nPh = nnz(model.getActivePhases());

                ds = zeros(numel(so), nPh);
                phIndices = model.getPhaseIndices();
                if model.water
                    ds(:, phIndices(1)) = dsw;
                end
                if model.oil
                    ds(:, phIndices(2)) = dso;
                end
                if model.gas
                    ds(:, phIndices(3)) = dsg;
                end
                if model.solvent
                    ds(:, phIndices(4)) = dss;
                end

                state = model.updateStateFromIncrement(state, ds, problem, 's', inf, model.dsMaxAbs);
                % We should *NOT* be solving for oil saturation for this to make sense
                assert(~any(strcmpi(vars, 'so')));

                state = computeFlashBlackOilSolvent(state, state0, model, st);
                state.s = bsxfun(@rdivide, state.s, sum(state.s, 2));

                %  We have explicitly dealt with rs/rv properties, remove from list
                %  meant for autoupdate.
                [vars, ix] = model.stripVars(vars, {'sw', 'so', 'sg', 'ss', 'rs', 'rv', 'x'});
                removed(~removed) = removed(~removed) | ix;

            end

            % We may have solved for a bunch of variables already if we had
            % disgas / vapoil enabled, so we remove these from the
            % increment and the linearized problem before passing them onto
            % the generic reservoir update function.
            problem.primaryVariables = vars;
            dx(removed) = [];

            % Parent class handles almost everything for us
            [state, report] = updateState@ReservoirModel(model, state, problem, dx, drivingForces);        
        end

        % --------------------------------------------------------------------%
        function [isActive, phInd] = getActivePhases(model)
            % Get active flag for the canonical phase ordering (water, oil
            % gas as on/off flags).
            isActive = [model.water, model.oil, model.gas, model.solvent];
            if nargout > 1
                phInd = find(isActive);
            end
        end

        % --------------------------------------------------------------------%
        function state = storeFluxes(model, state, vW, vO, vG, vS)
            % Utility function for storing the interface fluxes in the state
            isActive = model.getActivePhases();

            internal = model.operators.internalConn;
            state.flux = zeros(numel(internal), sum(isActive));
            phasefluxes = {value(vW), value(vO), value(vG), value(vS)};
            state = model.setPhaseData(state, phasefluxes, 'flux', internal);
        end

         % --------------------------------------------------------------------%
        function state = storeMobilities(model, state, mobW, mobO, mobG, mobS)
            % Utility function for storing the mobilities in the state
            isActive = model.getActivePhases();

            state.mob = zeros(model.G.cells.num, sum(isActive));
            mob = {value(mobW), value(mobO), value(mobG), value(mobS)};
            state = model.setPhaseData(state, mob, 'mob');
        end

        % --------------------------------------------------------------------%
        function state = storeUpstreamIndices(model, state, upcw, upco, upcg, upcs)
            % Store upstream indices, so that they can be reused for other
            % purposes.
            isActive = model.getActivePhases();

            nInterfaces = size(model.operators.N, 1);
            state.upstreamFlag = false(nInterfaces, sum(isActive));
            mob = {upcw, upco, upcg, upcs};
            state = model.setPhaseData(state, mob, 'upstreamFlag');
        end

        % --------------------------------------------------------------------%
        function state = storeDensity(model, state, rhoW, rhoO, rhoG, rhoS)
            % Store compressibility / surface factors for plotting and
            % output.
            isActive = model.getActivePhases();

            state.rho = zeros(model.G.cells.num, sum(isActive));
            rho = {value(rhoW), value(rhoO), value(rhoG), value(rhoS)};
            state = model.setPhaseData(state, rho, 'rho');
        end

        % --------------------------------------------------------------------%
        function state = storebfactors(model, state, bW, bO, bG, bS)
            % Store compressibility / surface factors for plotting and
            % output.
            isActive = model.getActivePhases();

            state.bfactor = zeros(model.G.cells.num, sum(isActive));
            b = {value(bW), value(bO), value(bG), value(bS)};
            state = model.setPhaseData(state, b, 'bfactor');
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