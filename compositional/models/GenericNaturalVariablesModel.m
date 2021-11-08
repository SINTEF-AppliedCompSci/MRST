classdef GenericNaturalVariablesModel < NaturalVariablesCompositionalModel & GenericReservoirModel
    properties
        
    end
    
    methods
        function model = GenericNaturalVariablesModel(varargin)
            model = model@NaturalVariablesCompositionalModel(varargin{:});
            model.OutputStateFunctions = {'ComponentTotalMass'};
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = getEquations@ReservoirModel(model, state0, state, dt, drivingForces, varargin{:});
        end
        
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces)
            % Discretize
            [eqs, flux, names, types] = model.FlowDiscretization.componentConservationEquations(model, state, state0, dt);
            src = model.FacilityModel.getComponentSources(state);
            [pressures, sat, mob, rho, X] = model.getProps(state, 'PhasePressures', 's', 'Mobility', 'Density', 'ComponentPhaseMassFractions');
            comps = cellfun(@(x, y) {x, y}, X(:, model.getLiquidIndex), X(:, model.getVaporIndex), 'UniformOutput', false);
            
            
            eqs = model.addBoundaryConditionsAndSources(eqs, names, types, state, ...
                                                             pressures, sat, mob, rho, ...
                                                             {}, comps, ...
                                                             drivingForces);
            
            % Add sources
            eqs = model.insertSources(eqs, src);
            % Assemble equations
            for i = 1:numel(eqs)
                eqs{i} = model.operators.AccDiv(eqs{i}, flux{i});
            end
            % Natural variables part
            n_hc = model.EOSModel.getNumberOfComponents();
            twoPhase = model.getTwoPhaseFlag(state);
            cnames = model.EOSModel.getComponentNames;
            f = model.getProps(state, 'Fugacity');
            f_eqs = cell(1, n_hc);
            f_names = cell(1, n_hc);
            f_types = cell(1, n_hc);
            
            
            for i = 1:n_hc
                if any(twoPhase)
                    f_eqs{i} = (f{i, 1}(twoPhase) - f{i, 2}(twoPhase))/barsa;
                end
                f_names{i} = ['f_', cnames{i}];
                f_types{i} = 'fugacity';
            end
            % Get facility equations
            [weqs, wnames, wtypes, state] = model.FacilityModel.getModelEquations(state0, state, dt, drivingForces);
            % Finally assemble
            eqs = [eqs, weqs, f_eqs];
            names = [names, wnames, f_names];
            types = [types, wtypes, f_types];
        end
        
        function names = getComponentNames(model)
            names = cellfun(@(x) x.name, model.Components, 'UniformOutput', false);
        end

        function [state, report] = updateState(model, state, problem, dx, forces)
            [state, report] = updateState@NaturalVariablesCompositionalModel(model, state, problem, dx, forces);
            if ~isempty(model.FacilityModel)
                state = model.FacilityModel.applyWellLimits(state);
            end
        end
        
        function model = validateModel(model, varargin)
            % Validate model.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.validateModel`
            if isempty(model.FacilityModel) || ~isa(model.FacilityModel, 'GenericFacilityModel')
                model.FacilityModel = GenericFacilityModel(model);
            end
            if isempty(model.Components)
                names_eos = model.EOSModel.getComponentNames();
                % Add in additional immiscible phases
                [sn_regular, phases_regular] = model.getNonEoSPhaseNames();
                nreg = numel(sn_regular);
                enames = cell(1, nreg);
                for i = 1:nreg
                    enames{i} = phases_regular{i};
                end
                names = [names_eos, enames];
                
                nc = numel(names);
                model.Components = cell(1, nc);
                p = model.FacilityModel.pressure;
                T = model.FacilityModel.T;
                for ci = 1:nc
                    name = names{ci};
                    switch name
                        case {'water', 'oil', 'gas'}
                            ix = model.getPhaseIndex(upper(name(1)));
                            c = ImmiscibleComponent(name, ix);
                        otherwise
                            c = getEOSComponent(model, p, T, name, ci);
                    end
                    model.Components{ci} = c;
                end
            end
            model = validateModel@NaturalVariablesCompositionalModel(model, varargin{:});
        end
        
        function model = setupStateFunctionGroupings(model, varargin)
            model = setupStateFunctionGroupings@NaturalVariablesCompositionalModel(model, varargin{:});
            model.FlowDiscretization.GravityPotentialDifference.saturationWeighting = true;
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@ThreePhaseCompositionalModel(model, state0, state, dt, drivingForces);
            if model.outputFluxes
                f = model.getProp(state, 'PhaseFlux');
                nph = numel(f);
                state.flux = zeros(model.G.faces.num, nph);
                state.flux(model.operators.internalConn, :) = value(f);
                if ~isempty(drivingForces.bc)
                    [p, s, mob, rho, b] = model.getProps(state, 'PhasePressures', 's', 'Mobility', 'Density', 'ShrinkageFactors');
                    sat = expandMatrixToCell(s);
                    [~, ~, ~, fRes] = getBoundaryConditionFluxesAD(model, p, sat, mob, rho, b, drivingForces.bc);
                    idx = model.getActivePhases();
                    fWOG = cell(3, 1);
                    fWOG(idx) = fRes;

                    state = model.storeBoundaryFluxes(state, fWOG{1}, fWOG{2}, fWOG{3}, drivingForces);
                end
            end
        end

        function [vars, names, origin] = getPrimaryVariables(model, state)
            % Get primary variables from state, before a possible
            % initialization as AD.
            % Properties at current timestep
            [p, sL, sV, x, y] = model.getProps(state, ...
                'pressure', 'sL', 'sV', 'x', 'y');
            [pureLiquid, pureVapor, twoPhase] = model.getFlag(state);

            if 1
                stol = 1e-6;
                pureNonEoS = sL + sV < stol;
                sL(~pureVapor & pureNonEoS) = stol;
                sV(~pureLiquid & pureNonEoS) = stol;
            end
            z_tol = model.EOSModel.minimumComposition;

            x = ensureMinimumFraction(x, z_tol);
            y = ensureMinimumFraction(y, z_tol);
            x = expandMatrixToCell(x);
            y = expandMatrixToCell(y);

            ncomp = model.EOSModel.getNumberOfComponents();
            [xnames, ynames, cnames] = deal(model.EOSModel.getComponentNames());
            for i = 1:ncomp
                xnames{i} = ['v_', cnames{i}];
                ynames{i} = ['w_', cnames{i}];
            end
            n2ph = sum(twoPhase);
            
            if n2ph > 0
                wtmp = ones(n2ph, 1);
                w = cell(1, ncomp-1);
                [w{:}] = deal(wtmp);

                for i = 1:(ncomp-1)
                    w{i} = y{i}(twoPhase);
                end
                sl = sL(twoPhase);
            else
                w = cell(1, ncomp-1);
                sl = {[]};
            end

            if not(isempty(model.FacilityModel))
                [v, n, o] = model.FacilityModel.getPrimaryVariables(state);
            else
                [v, n, o] = deal({});
            end
            local_origin = class(model);
            
            component_names = xnames(1:end-1);
            comps = x(1:end-1);
            extra = model.getNonEoSPhaseNames();
            for i = 1:numel(extra)
                ns = ['s', extra(i)];
                component_names = [component_names, ns]; %#ok
                comps = [comps, model.getProp(state, ns)]; %#ok
            end
            vars = [p, comps, v, sl, w];
            names = ['pressure', component_names, n, 'sL', ynames(1:end-1)];

            offset = numel(component_names) + 1;
            origin = cell(1, numel(vars));
            [origin{:}] = deal(local_origin);
            origin((1:numel(v)) + offset) = o;
        end

        function state = initStateAD(model, state, vars, names, origin)
            [pureLiquid, pureVapor, twoPhase] = model.getFlag(state);
            singlePhase = pureLiquid | pureVapor;
            twoPhaseIx = find(twoPhase);
            nvar = numel(vars);
            removed = false(size(vars));
            cellJacMap = cell(nvar, 1);
            s = getSampleAD(vars{:});
            
            is_sl = strcmpi(names, 'sl');
            is_sv = strcmpi(names, 'sv');
                        
            cnames = model.EOSModel.getComponentNames;
            ncomp = numel(cnames);
            x = cell(1, ncomp);
            w = cell(1, ncomp);
            y = cell(1, ncomp);
            x{end} = ones(model.G.cells.num, 1);
            w{end} = ones(sum(twoPhase), 1);
            
            for i = 1:ncomp-1
                name = cnames{i};
                is_x = strcmp(names, ['v_', name]);
                is_w = strcmp(names, ['w_', name]);
                x{i} = vars{is_x};
                x{end} = x{end}-x{i};
                if any(twoPhase)
                    w{i} = vars{is_w};
                    w{end} = w{end}-w{i};
                end
                removed(is_x | is_w) = true;
                % We note that this is a two-phase subset
                cellJacMap{is_w} = twoPhaseIx;
            end

            for i = 1:ncomp
                y{i} = singlePhase.*x{i};
                if any(twoPhase)
                    y{i}(twoPhase) = w{i};
                end
            end
            state = model.setProps(state, ...
                {'liquidMoleFractions', 'vaporMoleFractions'}, {x, y});
            % Deal with saturations
            phases = model.getPhaseNames();
            nph = numel(phases);
            sE = zeros(model.G.cells.num, 1);
            sat = cell(1, nph);
            if nph > 2
                extra = model.getNonEoSPhaseNames();
                ne = numel(extra);
                for i = 1:ne
                    e = extra(i);
                    is_s = strcmpi(names, ['s', e]);
                    S = vars{is_s};
                    % Add to sat
                    sE = sE + S;
                    sat{phases == e} = S;
                    removed(is_s) = true;
                end
            end
            % EoS values
            [sL, sV] = model.getProps(state, 'sL', 'sV');
            sL = model.AutoDiffBackend.convertToAD(sL, s);
            sV = model.AutoDiffBackend.convertToAD(sV, s);
            if any(twoPhase)
                if nph > 2
                    if any(pureLiquid)
                        sL(pureLiquid) = 1 - sE(pureLiquid);
                    end
                    if any(pureVapor)
                        sV(pureVapor) = 1 - sE(pureVapor);
                    end
                end
                % Set liquid saturation in two-phase cells
                sl = vars{is_sl};
                sL(twoPhase) = sl;
                % Set vapor saturation in two-phase cells
                sV(twoPhase) = 1 - sl - sE(twoPhase);
                if nph == 2
                    cellJacMap{is_sl} = twoPhaseIx;
                end
            end
            removed(is_sv | is_sl) = true;

            if nph > 2
                [sL, sV] = model.setMinimumTwoPhaseSaturations(state, sE, sL, sV, pureLiquid, pureVapor);
            end
            sat{model.getLiquidIndex} = sL;
            sat{model.getVaporIndex} = sV;

            state = model.setProp(state, 's', sat);

            if ~isempty(model.FacilityModel)
                % Select facility model variables and pass them off to attached
                % class.
                fm = class(model.FacilityModel);
                isF = strcmp(origin, fm);
                state = model.FacilityModel.initStateAD(state, vars(isF), names(isF), origin(isF));
                removed = removed | isF;
            end
            state.cellJacMap = cellJacMap;
            % Set up state with remaining variables
            state = initStateAD@ReservoirModel(model, state, vars(~removed), names(~removed), origin(~removed));
        end

        function forces = validateDrivingForces(model, forces, varargin)
            forces = validateDrivingForces@NaturalVariablesCompositionalModel(model, forces, varargin{:});
            forces = validateCompositionalForces(model, forces, varargin{:});
        end
        
        function problem = setupLinearizedProblem(model, eqs, types, names, primaryVars, state, dt)
            s = eqs{1};
            if model.reduceLinearSystem && isa(s, 'ADI')
                problem = ReducedLinearizedSystem(eqs, types, names, primaryVars, state, dt);
                [~, ~, twoPhase] = model.getFlag(state);
                % Switch first composition and first saturation in
                % two-phase region to ensure invertible
                % Schur-complement
                twoPhaseIx = find(twoPhase);
                cn = model.EOSModel.getComponentNames();
                xInd = strcmpi(primaryVars, ['v_', cn{1}]);
                sInd = strcmpi(primaryVars, 'sL');
                offsets = cumsum([0; s.getNumVars()]);
                reorder = 1:offsets(end);
                start = offsets(xInd) + twoPhaseIx;
                stop = offsets(sInd) + (1:sum(twoPhase));
                reorder(start) = stop;
                reorder(stop) = start;
                problem.reorder = reorder;
                problem.keepNum = offsets(sInd);
            else
                problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
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
