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
    %   compFluid - CompositionalFluid instance describing the species
    %               present.
    %
    % RETURNS:
    %  model - Initialized class instance
    %
    % SEE ALSO:
    %   `ThreePhaseCompositionalModel`, `OverallCompositionCompositionalModel`
    properties
        allowLargeSaturations = false; % Allow sum of saturations larger than unity (experimental option)
        reduceLinearSystem = true; % Return ReducedLinearizedSystem instead of LinearizedSystemAD
        maxPhaseChangesNonLinear = inf; % Maximum number of phase transitions for a given cell, during a nonlinear step (experimental option)
    end
    
    methods
        function model = NaturalVariablesCompositionalModel(G, rock, fluid, compFluid, varargin)
            model = model@ThreePhaseCompositionalModel(G, rock, fluid, compFluid, varargin{:});
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsNaturalVariables(state0, state, model, dt, ...
                            drivingForces, varargin{:});
        end
        
        function [ds, dx, vars] = getSaturationIncrements(model, dx, vars, state)
            [pureLiquid, pureVapor, twoPhase] = model.getFlag(state);
            % ----------------------- GAS PHASE --------------------------%
            [dsG, dx, vars] = getSatUpdateInternal(model, 'satg', dx, vars, twoPhase);
            [dsO, dx, vars] = getSatUpdateInternal(model, 'sato', dx, vars, twoPhase);
            [dsW, dx, vars] = getSatUpdateInternal(model, 'satw', dx, vars, twoPhase);
            if model.water
                if ~any(strcmpi(vars, 'sGsO'))
                    dsO(pureLiquid) = -dsW(pureLiquid);
                    dsG(pureVapor) = -dsW(pureVapor);
                end
                ds = [dsW, dsO, dsG];
            else
                ds = [dsO, dsG];
            end
        end

        function [dx, dy, increments, vars] = getPhaseCompositionIncrements(model, increments, vars, state)
            [pureLiquid, pureVapor, twoPhase] = model.getFlag(state);
            
            ncomp = model.EOSModel.fluid.getNumberOfComponents();
            cnames = model.EOSModel.fluid.names;
            
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
            
            [pureLiquid, pureVapor, twoPhase] = model.getFlag(state);
            mixvar = strcmpi(vars, 'sGsO');
            if any(mixvar)
                % Sequential implicit stuff
                dMix = dx{mixvar};
                ds(~pureLiquid, end) = dMix(~pureLiquid);
                ds(pureLiquid, end-1) = dMix(pureLiquid);
                [vars, removed] = model.stripVars(vars, 'sGsO');
                dx = dx(~removed);
            end
            dL = [];
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            state0 = state;
            dx0 = dx;
            vars = problem.primaryVariables;
            if model.allowLargeSaturations
                dsMax = max(sum(state.s, 2), 1).*model.dsMaxAbs;
            else
                dsMax = model.dsMaxAbs;
            end
            [dxMax, dyMax] = deal(model.dzMaxAbs);
            % Get saturation and composition updates
            
            [ds, deltax, deltay, dL, dx, vars] = model.getHyperbolicUpdates(problem, dx, vars, state);
            % dsTol = 2*dsMax;
            % for i = 1:size(ds, 2)
            %    subs = abs(ds(:, i)) > dsTol;
            %    ds(subs, i) = dsTol(subs).*sign(ds(subs, i));
            % end
            w_sat = min(dsMax./max(abs(ds), [], 2), 1);
            w_x = min(dxMax./max(abs(deltax), [], 2), 1);
            w_y = min(dyMax./max(abs(deltay), [], 2), 1);
            w_xy = min(w_x, w_y);
            
            % Minimum relax factor
            w = min(w_sat, w_xy);
            % Update sat
            capunit = @(x) min(max(x, 0), 1);
            if any(ds(:) ~= 0)
                ds_relax = bsxfun(@times, w, ds);
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

            if model.water
                sT = sum(state.s, 2);
                space = (sT - state.s(:, 1))./sT;
            else
                space = 1;
            end
            state.dz = max(max(abs(bsxfun(@times, deltax, space)), bsxfun(@times, deltay, space)));
            dz2 = max(abs(bsxfun(@times, state.components - state0.components, space)));
            state.dz = max(state.dz, dz2);
        end

        function state = flashPhases(model, state, state0, s_uncap, xyUpdated, iteration)
            isTransport = isa(model, 'TransportNaturalVariablesModel');
            isPressure = isa(model, 'PressureNaturalVariablesModel');
            
            oilIndex = 1 + model.water;
            gasIndex = 2 + model.water;

            p = state.pressure;
            T = state.T;
            [isLiquid0, isVapor0, isTwoPh0] = model.getFlag(state);
            if isPressure
                x = state.x;
                y = state.y;
                z = bsxfun(@times, x, state.L) + bsxfun(@times, 1-state.L, y);
            else
                x = state.x;
                y = state.y;
                z = bsxfun(@times, x, state.L) + bsxfun(@times, 1-state.L, y);
                z(isLiquid0, :) = state.x(isLiquid0, :);
                z(isVapor0, :) = state.y(isVapor0, :);
                z = bsxfun(@rdivide, z, sum(z, 2));
            end
            
            if iteration == 1
                state.switchCount = 0*state.pressure;
            end
            
            if 0
                [stable, x, y] = phaseStabilityTest(model.EOSModel, z, p, T);
            else
                x = state.x;
                y = state.y;
                act = ~isTwoPh0;
                stable = act;
                [stable(act), x(act, :), y(act, :)] =...
                    model.EOSModel.performPhaseStabilityTest(p(act, :), T(act), z(act, :));
            end
            
            locked = state.switchCount > model.maxPhaseChangesNonLinear;
            stable(locked) = true;
            
            sO = state.s(:, oilIndex);
            sG = state.s(:, gasIndex);
            tol = 1e-8;

            toEpsOil = isVapor0 & ~stable;
            toEpsGas = isLiquid0 & ~stable;
            
            sO_uncap = s_uncap(:, oilIndex);
            sG_uncap = s_uncap(:, gasIndex);
            if model.water
                sW = state.s(:, 1);
            else
                sW = 0*sO;
            end
            toOnlyOil = isTwoPh0 & sG_uncap <= 0;
            toOnlyGas = isTwoPh0 & sO_uncap <= 0;

            bad = toOnlyOil & toOnlyGas;
            if any(bad)
                gasLargest = abs(state0.s(:, oilIndex)) < abs(state0.s(:, gasIndex));
                toOnlyGas(bad &  gasLargest) = true;
                toOnlyGas(bad & ~gasLargest) = false;
                toOnlyOil(bad & ~gasLargest) = true;
                toOnlyOil(bad &  gasLargest) = false;
            end
            
            isPureLiquid = (stable & isLiquid0) | toOnlyOil;
            isPureVapor  = (stable & isVapor0 & ~isPureLiquid)  | toOnlyGas;

            switched = toEpsOil | toEpsGas;
            if any(switched)
                state.x(switched, :) = x(switched, :);
                state.y(switched, :) = y(switched, :);
            end
            state.switchCount = state.switchCount + double(switched | toOnlyGas | toOnlyOil);
            if isa(model, 'TransportNaturalVariablesModel')
                sMax = sum(state.s, 2);
            else
                sMax = ones(model.G.cells.num, 1);
            end
            if model.water
                sMax = sMax - sW;
            end
            
            ds_oswitch = tol.*sMax(toEpsOil);
            ds_gswitch = tol.*sMax(toEpsGas);
            
            sG(toEpsOil) = sMax(toEpsOil) - ds_oswitch;
            sO(toEpsGas) = sMax(toEpsGas) - ds_gswitch;

            sO(toEpsOil) = ds_oswitch;
            sG(toEpsGas) = ds_gswitch;
            

            sO(isPureVapor) = 0;
            sG(isPureLiquid) = 0;
            if isPressure
                sO(isPureLiquid) = sMax(isPureLiquid);
                sG(isPureVapor) = sMax(isPureVapor);
            else
                sO(toOnlyOil) = sMax(toOnlyOil);
                sG(toOnlyGas) = sMax(toOnlyGas);
            end

            if model.water
                state.s = [sW, sO, sG];
            else
                state.s = [sO, sG];
            end
            state = model.setFlag(state, isPureLiquid, isPureVapor);
            
            if isPressure
                pure = isPureVapor | isPureLiquid;
                state.x(pure, :) = state.z0(pure, :);
                state.y(pure, :) = state.z0(pure, :);
            else
                state.x(isPureVapor, :) = state.y(isPureVapor, :);
                state.y(isPureLiquid, :) = state.x(isPureLiquid, :);                
            end

            state = model.updateSecondaryProperties(state);
            
            [isLiquid, isVapor, isTwoPh] = model.getFlag(state);
            state.switched = (isTwoPh0 ~= isTwoPh);

            state.x = ensureMinimumFraction(state.x);
            state.y = ensureMinimumFraction(state.y);
            if any(xyUpdated)
                state.components = bsxfun(@times, state.x, state.L) + bsxfun(@times, state.y, (1-state.L));
            end
            bad = all(state.s == 0, 2);
            if any(bad)
                unset = sMax;
                state.s(bad & isPureVapor, gasIndex) = unset(bad & isPureVapor);
                state.s(bad & isPureLiquid, oilIndex) = unset(bad & isPureLiquid);
                state.s(bad & isTwoPh, gasIndex) = unset(bad & isTwoPh)/2;
                state.s(bad & isTwoPh, oilIndex) = unset(bad & isTwoPh)/2;
            end
            dispif(model.verbose > 1, '%d 2ph, %d oil, %d gas [%d switched, %d locked]\n', ...
                nnz(state.flag == 0), nnz(state.flag == 1), nnz(state.flag == 2), nnz(state.switched), nnz(locked));
            
        end
        
        function state = updateSecondaryProperties(model, state)
            sO = state.s(:, 1 + model.water);
            sG = state.s(:, 2 + model.water);
            sO = sO./(sO + sG);
            sG = sG./(sO + sG);
            x = state.x;
            y = state.y;
            z = state.components;
            p = state.pressure;
            temp = state.T;

            eos = model.EOSModel;
            [Z_L, Z_V] = eos.getProperties(p, temp, x, y, z, sO, sG);
            if isa(model, 'PressureNaturalVariablesModel')
                twoPhase = state.flag == 0;
                [~, mv] = max(state.z0, [], 2);
                
                mv = sub2ind(size(state.z0), (1:numel(mv))', mv);
                L_new = (state.z0(mv) - state.y(mv))./(state.x(mv) - state.y(mv));
                state.L(twoPhase) = L_new(twoPhase);
            else
                rhoO = model.PropertyModel.computeMolarDensity(p, x, Z_L, temp, true);
                rhoG = model.PropertyModel.computeMolarDensity(p, y, Z_V, temp, false);
                L = rhoO.*sO./(rhoO.*sO + rhoG.*sG);
                state.L = double(L);
            end
            
            state.Z_L = Z_L;
            state.Z_V = Z_V;
        end

        
        function [xM,  yM,  rhoO,  rhoG,  muO,  muG, f_L, f_V,...
                  xM0, yM0, rhoO0, rhoG0] = ...
                  getTimestepPropertiesEoS(model, state, state0, p, temp, x, y, z, sO, sG, cellJacMap)
            eos = model.EOSModel;
            [Z_L, Z_V, f_L, f_V] = eos.getProperties(p, temp, x, y, z, sO, sG, state, cellJacMap);

            [xM, rhoO, muO] = model.getFlowPropsNatural(p, x, Z_L, temp, true);
            [yM, rhoG, muG] = model.getFlowPropsNatural(p, y, Z_V, temp, false);
            
            [p0, temp0, x0, y0] = model.getProps(state0, 'pressure', 'T', 'x', 'y');
            x0 = expandMatrixToCell(x0);
            y0 = expandMatrixToCell(y0);

            [xM0, rhoO0] = model.getFlowPropsNatural(p0, x0, state0.Z_L, temp0, true);
            [yM0, rhoG0] = model.getFlowPropsNatural(p0, y0, state0.Z_V, temp0, false);

        end
        
        function [xM, rho, mu] = getFlowPropsNatural(model, p, x, Z, T, isLiquid)
            xM = model.EOSModel.getMassFraction(x);
            rho = model.PropertyModel.computeDensity(p, x, Z, T, isLiquid);
            if nargout > 2
                mu = model.PropertyModel.computeViscosity(p, x, Z, T, isLiquid);
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
    end
    
    methods(Access=protected)
       function [ds, dx, vars] = getSatUpdateInternal(model, name, dx, vars, twoPhase)
        dsgix = strcmpi(vars, name);
        if any(dsgix)
            ds_tmp = dx{dsgix};
            if numel(ds_tmp) == nnz(twoPhase)
                ds = zeros(model.G.cells.num, 1);
                ds(twoPhase) = ds_tmp;
            else
                assert(numel(ds_tmp) == model.G.cells.num);
                ds = ds_tmp;
            end
            [vars, removed] = model.stripVars(vars, name);
            dx = dx(~removed);
        else
            ds = zeros(model.G.cells.num, 1);
        end
    end
 
    end
end


