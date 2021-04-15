classdef NaturalVariablesCompositionalModel < ThreePhaseCompositionalModel
    % Natural variables model for compositional problems
    %
    % SYNOPSIS:
    %   model = NaturalVariablesCompositionalModel(G, rock, fluid, compFluid)
    %
    % DESCRIPTION:
    %   The natural variables model relies on separate primary variables
    %   for saturations in the multiphase region. This makes it possible to
    %   perform certain heuristics to ensure convergence (e.g. saturation
    %   chopping) since there is a weaker correspondence between
    %   saturations and the compositions. However, this formulation
    %   traverses the phase boundary as a number of discrete steps and can
    %   take more nonlinear iterations for some problems.
    %
    % PARAMETERS:
    %   G         - Grid structure
    %   rock      - Rock structure for the reservoir
    %   fluid     - The flow fluid, containing relative permeabilities,
    %               surface densities and flow properties for the
    %               aqueous/water phase (if present)
    %   mixture   - CompositionalMixture instance describing the species
    %               present.
    %
    % RETURNS:
    %  model - Initialized class instance
    %
    % SEE ALSO:
    %   `ThreePhaseCompositionalModel`, `OverallCompositionCompositionalModel`
    properties
        allowLargeSaturations = false;  % Allow sum of saturations larger than unity (experimental option)
        reduceLinearSystem = true;      % Return ReducedLinearizedSystem instead of LinearizedSystemAD
        maxPhaseChangesNonLinear = inf; % Maximum number of phase transitions for a given cell, during a nonlinear step (experimental option)
        checkStableTransition = false;  % Do not transition to single-phase if stable
        saturationEpsilon = 1e-6;       % Epsilon for saturation during phase change
        flashFromSinglePhase = false;   % Perform flash for cells existing single-phase
    end
    
    methods
        function model = NaturalVariablesCompositionalModel(G, rock, fluid, mixture, varargin)
            model = model@ThreePhaseCompositionalModel(G, rock, fluid, mixture, varargin{:});
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsNaturalVariables(state0, state, model, dt, ...
                            drivingForces, varargin{:});
        end
        
        function [ds, dx, vars] = getSaturationIncrements(model, dx, vars, state)
            [pureLiquid, pureVapor, twoPhase] = model.getFlag(state);
            phases = model.getPhaseNames();
            nph = numel(phases);
            ds = zeros(model.G.cells.num, nph);
            isFound = false(1, nph);
            for i = 1:nph
                [ds(:, i), dx, vars, isFound(i)] = getSatUpdateInternal(model, phases(i), dx, vars, twoPhase);
            end
            missing = ~isFound;
            if sum(missing) == 1
                ds(:, missing) = -sum(ds(:, isFound), 2);
            end
            if nph > 2
                li = model.getLiquidIndex();
                vi = model.getVaporIndex();
                tmp = value(ds);
                ds(pureLiquid, li) = -sum(tmp(pureLiquid, :), 2);
                ds(pureVapor, vi) = -sum(tmp(pureVapor, :), 2);
            end
        end

        function [dx, dy, increments, vars] = getPhaseCompositionIncrements(model, increments, vars, state)
            [pureLiquid, pureVapor, twoPhase] = model.getFlag(state);
            
            ncomp = model.EOSModel.getNumberOfComponents();
            cnames = model.EOSModel.getComponentNames();
            
            [dx, dy] = deal(zeros(model.G.cells.num, ncomp));
            found = false(ncomp, 1);

            for i = 1:ncomp
                namev = ['v_', cnames{i}];
                namew = ['w_', cnames{i}];
                vix = strcmpi(vars, namev);
                wix = strcmpi(vars, namew);
                if any(wix)
                    found(i) = true;
                    dx(:, i) = ~pureVapor.*increments{vix};
                    dy(pureVapor, i) = increments{vix}(pureVapor);
                    dy(twoPhase, i) = increments{wix};
                    [vars, removed] = model.stripVars(vars, {namev, namew});
                    increments = increments(~removed);
                end
            end
            
            n_found = nnz(found);
            if n_found == ncomp-1
                dx(:, ~found) = -sum(dx(:, found), 2);
                dy(:, ~found) = -sum(dy(:, found), 2);
            else
                assert(n_found == 0 || n_found == ncomp);
            end
        end
        
        function [ds, deltax, deltay, dL, dx, vars] = getHyperbolicUpdates(model, problem, dx, vars, state)
            [ds, dx, vars] = model.getSaturationIncrements(dx, vars, state);
            [deltax, deltay, dx, vars] = model.getPhaseCompositionIncrements(dx, vars, state);
            
            pureLiquid = model.getFlag(state);
            mixvar = strcmpi(vars, 'sGsO');
            if any(mixvar)
                % Sequential implicit stuff
                dMix = dx{mixvar};
                ds(~pureLiquid, model.getVaporIndex()) = dMix(~pureLiquid);
                ds(pureLiquid, model.getLiquidIndex()) = dMix(pureLiquid);
                [vars, removed] = model.stripVars(vars, 'sGsO');
                dx = dx(~removed);
            end
            dL = [];
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            state0 = state;
            vars = problem.primaryVariables;
            if model.allowLargeSaturations
                dsMax = max(sum(state.s, 2), 1).*model.dsMaxAbs;
            else
                dsMax = model.dsMaxAbs;
            end
            [dxMax, dyMax] = deal(model.dzMaxAbs);
            % Get saturation and composition updates
            
            [ds, deltax, deltay, dL, dx, vars] = model.getHyperbolicUpdates(problem, dx, vars, state);
            w_sat = min(dsMax./max(abs(ds), [], 2), 1);
            w_x = min(dxMax./max(abs(deltax), [], 2), 1);
            w_y = min(dyMax./max(abs(deltay), [], 2), 1);
            w_xy = min(w_x, w_y);
            
            % Minimum relax factor
            w = min(w_sat, w_xy);
            % Update sat
            capunit = @(x) min(max(x, 0), 1);
            nph = model.getNumberOfPhases();
            isEosPhase = true(1, nph);
            isEoSPhase(model.getEoSPhaseIndices) = true;

            if any(ds(:) ~= 0)
                ds_relax = ds;
                ds_relax(:, ~isEoSPhase) = ds(:, ~isEoSPhase).*min(dsMax./max(abs(ds(:, ~isEoSPhase)), [], 2), 1);
                ds_relax(:, isEosPhase) = bsxfun(@times, ds(:, isEosPhase), w);
                s_uncap = state.s + ds_relax;
                if model.allowLargeSaturations
                    state.s = max(s_uncap, 0);
                else
                    state.s = capunit(s_uncap);
                    state.s = bsxfun(@rdivide, state.s, sum(state.s, 2));
                end
            else
                s_uncap = state.s;
            end

            % Update dx, dy
            if any(deltax(:) ~= 0) || any(deltay(:) ~= 0)
                state.x = max(state.x + bsxfun(@times, deltax, w), 0);
                state.y = max(state.y + bsxfun(@times, deltay, w), 0);
                
                state.x = capunit(state.x);
                state.x = bsxfun(@rdivide, state.x, sum(state.x, 2));

                state.y = capunit(state.y);
                state.y = bsxfun(@rdivide, state.y, sum(state.y, 2));
                % Update K-values
                isTwoPhase = model.getTwoPhaseFlag(state);
                state.x = ensureMinimumFraction(state.x, model.EOSModel.minimumComposition);
                state.y = ensureMinimumFraction(state.y, model.EOSModel.minimumComposition);

                state.K(isTwoPhase, :) = state.y(isTwoPhase, :)./state.x(isTwoPhase, :);
                xyUpdated = true;
            else
                xyUpdated = false;
            end
            
            if ~isempty(dL)
                state.L = capunit(state.L + w.*dL);
            end
            
            problem.primaryVariables = vars;
            [state, report] = updateState@ThreePhaseCompositionalModel(model, state, problem, dx, drivingForces);
            state = model.flashPhases(state, state0, s_uncap, xyUpdated, problem.iterationNo);
            sT = sum(state.s, 2);
            sT_hc = sum(state.s(:, isEoSPhase), 2);
            s_hc = sT_hc./sT;
            dx = model.computeChange(deltax, s_hc);
            dy = model.computeChange(deltay, s_hc);
            dz = model.computeChange(state.components - state0.components, s_hc);
            state.dz = max(dz, max(dx, dy));
        end

        function state = flashPhases(model, state, state0, s_uncap, xyUpdated, iteration)
            isTransport = isa(model, 'TransportNaturalVariablesModel');
            isPressure = isa(model, 'PressureNaturalVariablesModel');
            
            liquidIndex = model.getLiquidIndex();
            vaporIndex = model.getVaporIndex();
            % Get pressure, temperature and compositions
            [p, T, K, x, y] = model.getProps(state, 'pressure', 'temperature', 'K', 'x', 'y');
            % Get indicators for phase state at linearization which gave us
            % updates. The pure liquid or pure vapor cells are candidates
            % for a transition to the two-phase region.
            [isLiquid0, isVapor0, isTwoPh0] = model.getFlag(state);
            % Compute z
            z = bsxfun(@times, x, state.L) + bsxfun(@times, 1-state.L, y);
            if ~isPressure
                z = bsxfun(@times, x, state.L) + bsxfun(@times, 1-state.L, y);
                z(isLiquid0, :) = state.x(isLiquid0, :);
                z(isVapor0, :) = state.y(isVapor0, :);
                z = bsxfun(@rdivide, z, sum(z, 2));
            end
            if iteration == 1
                state.switchCount = 0*state.pressure;
            end
            % Check phase stability for single-phase cells
            act = ~isTwoPh0;
            stable = act;
            [stable(act), x(act, :), y(act, :)] =...
                model.EOSModel.performPhaseStabilityTest(p(act, :), T(act, :), z(act, :), K(act, :));
            flashOut = model.flashFromSinglePhase;
            if flashOut
                fi = ~stable & act;
                substate = struct('pressure', p(fi), 'T', T(fi), 'components', z(fi, :), ...
                                  'K', K(fi, :), 's', state.s(fi, :), 'L', state.L(fi));
                substate = model.computeFlash(substate);
            end
            % Special check - we lock cells in single-phase region if a
            % slightly dangerous debug option is enabled.
            locked = state.switchCount > model.maxPhaseChangesNonLinear;
            stable(locked) = true;
            % Begun calculating new saturations
            sL = state.s(:, liquidIndex);
            sV = state.s(:, vaporIndex);
            tol = model.saturationEpsilon;
            % We introduce a small amount of the previously missing fluid
            % in the unstable cells
            toEpsLiq = isVapor0 & ~stable & ~flashOut;
            toEpsVap = isLiquid0 & ~stable & ~flashOut;
            % Un-capped saturation (indicating phase transition to
            % single-phase state)
            sL_uncap = s_uncap(:, liquidIndex);
            sV_uncap = s_uncap(:, vaporIndex);
            s_nonhc = 0;
            nhc = model.getNonEoSPhaseNames();
            for i = 1:numel(nhc)
                s_nonhc = s_nonhc + model.getProp(state, ['s', nhc(i)]);
            end
            % Transition to single-phase
            toOnlyLiq = isTwoPh0 & sV_uncap <= 0;
            toOnlyVap = isTwoPh0 & sL_uncap <= 0;
            % Optionally check if we are transitionining INTO a stable
            % state.
            if model.checkStableTransition
                toPure = toOnlyLiq | toOnlyVap;
                stable_next = false(size(toOnlyVap));
                stable_next(toPure) = model.EOSModel.performPhaseStabilityTest(p(toPure, :), T(toPure, :), z(toPure, :), K(toPure, :));
                
                badVap = ~stable_next & toOnlyVap;
                badLiq = ~stable_next & toOnlyLiq;
                
                toOnlyVap(badVap) = false;
                toOnlyLiq(badLiq) = false;
                
                toEpsLiq(badVap) = true;
                toEpsVap(badLiq) = true;
            end
            % Something very strange is going on - pick the liquid or vapor
            % state arbitrarily
            bad = toOnlyLiq & toOnlyVap;
            if any(bad)
                vapLargest = abs(state0.s(:, liquidIndex)) < abs(state0.s(:, vaporIndex));
                toOnlyVap(bad &  vapLargest) = true;
                toOnlyVap(bad & ~vapLargest) = false;
                toOnlyLiq(bad & ~vapLargest) = true;
                toOnlyLiq(bad &  vapLargest) = false;
            end
            % New phase flags
            isPureLiquid = (stable & isLiquid0) | toOnlyLiq;
            isPureVapor  = (stable & isVapor0 & ~isPureLiquid)  | toOnlyVap;
            % Cells switched to two-phase
            switched_to_twophase = toEpsLiq | toEpsVap;
            if any(switched_to_twophase)
                x_switched = x(switched_to_twophase, :);
                y_switched = y(switched_to_twophase, :);
                % Update K-value
                state.K(switched_to_twophase, :) = y_switched./x_switched;
                % Insert new values
                state.x(switched_to_twophase, :) = x_switched;
                state.y(switched_to_twophase, :) = y_switched;
            end
            state.switchCount = state.switchCount + double(switched_to_twophase | toOnlyVap | toOnlyLiq);
            % What is the maximum allowable total saturation?
            if model.allowLargeSaturations
                sMax = sum(state.s, 2);
            else
                sMax = ones(model.G.cells.num, 1);
            end
            % From this, define the minimum saturation
            sMin = tol.*sMax;
            % Ensure that there is a little bit of everything if
            % non-EoS phases exist
            sMax = sMax - s_nonhc;
            sMax = max(sMax, 1e-8);

            ds_oswitch = sMin(toEpsLiq);
            ds_gswitch = sMin(toEpsVap);
            
            sL(toEpsVap) = sMax(toEpsVap) - ds_gswitch;
            sV(toEpsLiq) = sMax(toEpsLiq) - ds_oswitch;

            sL(toEpsLiq) = ds_oswitch;
            sV(toEpsVap) = ds_gswitch;
            % Set single-phase saturations
            sL(isPureVapor) = 0;
            sV(isPureLiquid) = 0;
            if isTransport
                % Specific logic for transport model
                sL(toOnlyLiq) = sMax(toOnlyLiq);
                sV(toOnlyVap) = sMax(toOnlyVap);
            else
                sL(isPureLiquid) = sMax(isPureLiquid);
                sV(isPureVapor) = sMax(isPureVapor);
            end
            state = model.setProp(state, ['s', model.vaporPhase], sV);
            state = model.setProp(state, ['s', model.liquidPhase], sL);

            state = model.setFlag(state, isPureLiquid, isPureVapor);
            
            if isPressure
                % Ensure total mole fractions are kept fixed
                pure = isPureVapor | isPureLiquid;
                state.x(pure, :) = state.z0(pure, :);
                state.y(pure, :) = state.z0(pure, :);
            else
                % Going to single-phase implies that we end up with the
                % correct value.
                state.x(isPureVapor, :) = state.y(isPureVapor, :);
                state.y(isPureLiquid, :) = state.x(isPureLiquid, :);                
            end
            state = model.updateSecondaryProperties(state);
            isTwoPh = model.getTwoPhaseFlag(state);
            state.switched = (isTwoPh0 ~= isTwoPh);
            % Enforce minimum composition
            state.x = ensureMinimumFraction(state.x, model.EOSModel.minimumComposition);
            state.y = ensureMinimumFraction(state.y, model.EOSModel.minimumComposition);
            if any(xyUpdated)
                % Finally update total composition to be consistent with
                % liquid fraction
                state.components = bsxfun(@times, state.x, state.L) + bsxfun(@times, state.y, (1-state.L));
            end
            if flashOut
                fn = fieldnames(substate);
                for i = 1:numel(fn)
                    f = fn{i};
                    v = substate.(f);
                    if isnumeric(v)
                        state.(f)(fi, :) = v;
                    end
                end
            end
            % Check if we have any disasterous cells. Should not really
            % happen.
            bad = all(state.s == 0, 2);
            if any(bad)
                unset = sMax;
                state.s(bad & isPureVapor, vaporIndex) = unset(bad & isPureVapor);
                state.s(bad & isPureLiquid, liquidIndex) = unset(bad & isPureLiquid);
                state.s(bad & isTwoPh, vaporIndex) = unset(bad & isTwoPh)/2;
                state.s(bad & isTwoPh, liquidIndex) = unset(bad & isTwoPh)/2;
            end
            if model.verbose > 1
                fprintf('%d 2ph, %d oil, %d gas [%d switched, %d locked]\n', ...
                nnz(state.flag == 0), nnz(state.flag == 1), nnz(state.flag == 2), nnz(state.switched), nnz(locked));
            end
        end
        
        function state = updateSecondaryProperties(model, state)
            sL = state.s(:, model.getLiquidIndex());
            sV = state.s(:, model.getVaporIndex());
            sT = sL + sV;
            sL = sL./sT;
            sV = sV./sT;
            
            % Only non-EoS-phases
            missingPhases = sT == 0;
            sL(missingPhases) = state.L(missingPhases);
            sV(missingPhases) = 1-state.L(missingPhases);
            x = state.x;
            y = state.y;
            z = state.components;
            p = state.pressure;
            temp = state.T;

            eos = model.EOSModel;
            [Z_L, Z_V] = eos.getCompressibilityAndFugacity(p, temp, x, y, z, [], []);
            if isa(model, 'PressureNaturalVariablesModel')
                twoPhase = state.flag == 0;
                [~, mv] = max(state.z0, [], 2);
                
                mv = sub2ind(size(state.z0), (1:numel(mv))', mv);
                L_new = (state.z0(mv) - state.y(mv))./(state.x(mv) - state.y(mv));
                state.L(twoPhase) = L_new(twoPhase);
            else
                rhoO = model.PropertyModel.computeMolarDensity(eos, p, x, Z_L, temp, true);
                rhoG = model.PropertyModel.computeMolarDensity(eos, p, y, Z_V, temp, false);
                L = rhoO.*sL./(rhoO.*sL + rhoG.*sV);
                state.L = double(L);
            end
            
            state.Z_L = Z_L;
            state.Z_V = Z_V;
        end

        
        function [xM,  yM,  rhoL,  rhoV, muL,  muV, f_L, f_V,...
                  xM0, yM0, rhoL0, rhoV0] = ...
                  getTimestepPropertiesEoS(model, state, state0, p, temp, x, y, z, sO, sG, cellJacMap)
            eos = model.EOSModel;
            if isempty(eos.equilibriumConstantFunctions)
                [Z_L, Z_V, f_L, f_V] = eos.getCompressibilityAndFugacity(p, temp, x, y, z, state.Z_L, state.Z_V, cellJacMap);
            else
                [Z_L, Z_V] = eos.getCompressibilityAndFugacity(p, temp, x, y, z, state.Z_L, state.Z_V, cellJacMap);
                ncomp = numel(x);
                f_L = cell(1, ncomp);
                f_V = cell(1, ncomp);
                for i = 1:ncomp
                    K = eos.equilibriumConstantFunctions{i}(p, temp, z);
                    f_L{i} = K.*x{i};
                    f_V{i} = y{i};
                end
            end

            [xM, rhoL, muL] = model.getFlowPropsNatural(p, x, Z_L, temp, true);
            [yM, rhoV, muV] = model.getFlowPropsNatural(p, y, Z_V, temp, false);
            
            [p0, temp0, x0, y0] = model.getProps(state0, 'pressure', 'T', 'x', 'y');
            x0 = expandMatrixToCell(x0);
            y0 = expandMatrixToCell(y0);

            [xM0, rhoL0] = model.getFlowPropsNatural(p0, x0, state0.Z_L, temp0, true);
            [yM0, rhoV0] = model.getFlowPropsNatural(p0, y0, state0.Z_V, temp0, false);
        end
        
        function [xM, rho, mu] = getFlowPropsNatural(model, p, x, Z, T, isLiquid)
            eos = model.EOSModel;
            xM = eos.getMassFraction(x);
            rho = model.PropertyModel.computeDensity(eos, p, x, Z, T, isLiquid);
            if nargout > 2
                mu = model.PropertyModel.computeViscosity(eos, p, x, Z, T, isLiquid);
            end
        end
        
        function [convergence, values, names] = checkConvergence(model, problem, varargin)
            [convergence, values, names] = checkConvergence@ThreePhaseCompositionalModel(model, problem, varargin{:});
            
            if isfield(problem.state, 'switched') &&  ...
              (model.useIncTolComposition && ~isa(model, 'PressureNaturalVariablesModel'))
                % If we are using increment tolerances, we should also make
                % sure that no cells switched state from two-phase to
                % single-phase or vice versa
                convergence(1) = convergence(1) & ~any(problem.state.switched);
            end
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@ThreePhaseCompositionalModel(model, state0, state, dt, drivingForces);
            if isfield(state, 'mixing')
                state = rmfield(state, 'mixing');
            end
            if isfield(state, 'eos')
                state = rmfield(state, 'eos');
            end
        end
        
        function state = validateState(model, state)
            state = validateState@ThreePhaseCompositionalModel(model, state);
            [pureLiquid, pureVapor] = model.getFlag(state);
            singlePhase = pureLiquid | pureVapor;
            state.x(singlePhase, :) = state.components(singlePhase, :);
            state.y(singlePhase, :) = state.components(singlePhase, :);
        end
    end
    
    methods(Access=protected)
        function [ds, dx, vars, found] = getSatUpdateInternal(model, name, dx, vars, twoPhase)
            if strcmpi(name, model.liquidPhase)
                name = 'L';
            elseif strcmpi(name, model.vaporPhase)
                name = 'V';
            end
            sn = ['s', name];
            isdS = strcmpi(vars, sn);
            found = any(isdS);
            if found
                ds_tmp = dx{isdS};
                if numel(ds_tmp) == nnz(twoPhase)
                    ds = zeros(model.G.cells.num, 1);
                    ds(twoPhase) = ds_tmp;
                else
                    assert(numel(ds_tmp) == model.G.cells.num);
                    ds = ds_tmp;
                end
                [vars, removed] = model.stripVars(vars, sn);
                dx = dx(~removed);
                assert(any(removed))
            else
                ds = zeros(model.G.cells.num, 1);
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
