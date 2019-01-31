classdef ThreePhaseBlackOilModel < ReservoirModel
    % Three phase with optional dissolved gas and vaporized oil
properties
    % Flag deciding if gas can be dissolved into the oil phase
    disgas
    % Flag deciding if oil can be vaporized into the gas phase
    vapoil

    % Maximum relative Rs/Rv increment
    drsMaxRel
    % Maximum absolute Rs/Rv increment
    drsMaxAbs
end

methods
    function model = ThreePhaseBlackOilModel(G, rock, fluid, varargin)
        model = model@ReservoirModel(G, rock, fluid);

        % Typical black oil is disgas / dead oil, but all combinations
        % are supported
        model.vapoil = false;
        model.disgas = false;

        % Max increments
        model.drsMaxAbs = inf;
        model.drsMaxRel = inf;

        % Blackoil -> use CNV style convergence 
        model.useCNVConvergence = true;

        % All phases are present
        model.oil = true;
        model.gas = true;
        model.water = true;

        model = merge_options(model, varargin{:});

        d = model.inputdata;
        if ~isempty(d)
            % Assume ECL-style input deck, as this is the only
            % supported format at the moment.
            if isfield(d, 'RUNSPEC')
                if isfield(d.RUNSPEC, 'VAPOIL')
                    model.vapoil = d.RUNSPEC.VAPOIL;
                end
                if isfield(d.RUNSPEC, 'DISGAS')
                    model.disgas = d.RUNSPEC.DISGAS;
                end
            else
                error('Unknown dataset format!')
            end
        end
    end

    % --------------------------------------------------------------------%
    function model = validateModel(model, varargin)
        % Validate model.
        %
        % SEE ALSO:
        %   :meth:`ad_core.models.PhysicalModel.validateModel`

        if isempty(model.FlowPropertyFunctions)
            model.FlowPropertyFunctions = BlackOilFlowPropertyFunctions(model);
        end
        model = validateModel@ReservoirModel(model, varargin{:});
    end
    % --------------------------------------------------------------------%
    function [fn, index] = getVariableField(model, name, varargin)
        switch(lower(name))
            case {'rs', 'rv'}
                % RS and RV for gas dissolving into the oil phase and oil
                % components vaporizing into the gas phase respectively.
                fn = lower(name);
                index = 1;
            otherwise
                % Basic phases are known to the base class
                [fn, index] = getVariableField@ReservoirModel(model, name, varargin{:});
        end
    end
    
    function [vars, names, origin] = getPrimaryVariables(model, state)
        % Get primary variables from state
        % Properties at current timestep
        [p, sW, sG, rs, rv] = model.getProps(state, ...
            'pressure', 'water', 'gas', 'rs', 'rv');

        %Initialization of primary variables ----------------------------------
        st  = model.getCellStatusVO(state,  1-sW-sG, sW, sG);
        if model.disgas || model.vapoil
            % X is either Rs, Rv or Sg, depending on each cell's saturation status
            x = st{1}.*rs + st{2}.*rv + st{3}.*sG;
            gvar = 'x';
        else
            x = sG;
            gvar = 'sG';
        end
        vars = {p, sW, x};
        names = {'pressure', 'sW', gvar};
        origin = cell(1, numel(vars));
        [origin{:}] = deal(class(model));
        
        if not(isempty(model.FacilityModel))
            [v, n, o] = model.FacilityModel.getPrimaryVariables(state.wellSol);
            vars = [vars, v];
            names = [names, n];
            origin = [origin, o];
        end
    end

    function state = initStateAD(model, state, vars, names, origin)
        removed = false(size(vars));
        if model.disgas || model.vapoil
            isw = strcmpi(names, 'sw');
            isx = strcmpi(names, 'x');
            
            sW = vars{isw};
            x = vars{isx};
            sG = model.getProps(state, 'sg');
            st  = model.getCellStatusVO(state, 1-sW-sG, sW, sG);
            sG = st{2}.*(1-sW) + st{3}.*x;
            sO = st{1}.*(1-sW) + ~st{1}.*(1 - sW - sG);
            sat =  {sW, sO, sG};
            state = model.setProp(state, 's', sat);
            removed = removed | isx | isw;
        end
        if not(isempty(model.FacilityModel))
            % Temporary code...
            fm = class(model.FacilityModel);
            isF = strcmp(origin, fm);
            state = model.FacilityModel.initStateAD(state, vars(isF), names(isF), origin(isF));
%             state.FacilityState = struct('primaryVariables', {vars(isF)}, 'names', {names(isF)});
            removed = removed | isF;
        end
        
        % Set up state with remaining variables
        state = initStateAD@ReservoirModel(model, state, vars(~removed), names(~removed), origin(~removed));
        % Account for dissolution changing variables
        if model.disgas
            rsSat = model.getProp(state, 'RsMax');
            rs = ~st{1}.*rsSat + st{1}.*x;
            rs = rs.*(value(sO) > 0);
            state = model.setProp(state, 'rs', rs);
        end

        if model.vapoil
            rvSat = model.getProp(state, 'RvMax');
            rv = ~st{2}.*rvSat + st{2}.*x;
            rv = rv.*(value(sG) > 0);
            state = model.setProp(state, 'rv', rv);
        end
    end
    
    % --------------------------------------------------------------------%
    function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
        opt = struct('Verbose',     mrstVerbose,...
                    'reverseMode', false,...
                    'resOnly',     false,...
                    'iteration',   -1, ...
                    'drivingForces0', []);
        opt = merge_options(opt, varargin{:});
        
        
        % Define primary variables
        if opt.reverseMode
            [state0, primaryVars] = model.getReverseStateAD(state0);
            % The model must be validated with drivingForces so that the
            % FacilityModel gets updated.
            model = model.validateModel(drivingForces);
            state = model.getStateAD(state, false);
        else
            [state, primaryVars] = model.getStateAD(state, ~opt.resOnly);
        end

        
        [acc, divTerms, names, types] = conservationEquationsBlackOil(state0, state, model, dt, drivingForces);
        
        dissolved = model.getDissolutionMatrix(state.rs, state.rv);
        % Add in and setup well equations
        
        wellVars = state.FacilityState.primaryVariables;
        w = state.FacilityState.names;
        nw = numel(state.wellSol);
        wellMap = struct('isBHP', strcmp(w, 'bhp'), 'isRate', ...
                          ismember(w, {'qOs', 'qGs', 'qWs'}), ...
                         'extraMap', zeros(nw, nw));
        
        p = state.pressure;
        mob = state.FlowProps.Mobility;
        rho = state.FlowProps.Density;
        sat = state.s;
        pressures = state.FlowProps.PhasePressures;
        eqs = acc;

        [eqs, state] = model.addBoundaryConditionsAndSources(eqs, names, types, state, ...
                                                         pressures, sat, mob, rho, ...
                                                         dissolved, {}, ...
                                                         drivingForces);
        [eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, ...
                                                          names, types, ...
                                                          state0.wellSol, ...
                                                          state.wellSol, ...
                                                          wellVars, wellMap, ...
                                                          p, mob, rho, dissolved, ...
                                                          {}, dt, opt);
        for i = 1:numel(divTerms)
            eqs{i} = eqs{i} + divTerms{i};
        end
        state = value(state);
        problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
    end


    % --------------------------------------------------------------------%
    function state = validateState(model, state)
        % Check parent class
        state = validateState@ReservoirModel(model, state);
        nc = model.G.cells.num;
        if model.disgas
            % RS must be supplied for all cells. This may cause an error.
            model.checkProperty(state, 'rs', nc, 1);
            rsMax = model.getProp(state, 'rsMax');
            [sg, so, rs] = model.getProps(state, 'sg', 'so', 'rs');
            rs(sg > 0) = rsMax(sg > 0);
            rs = rs.*(so > 0);
            state = model.setProp(state, 'rs', rs);
        else
            % RS does not really matter. Assign single value.
            fn = model.getVariableField('rs');
            if ~isfield(state, fn)
                dispif(model.verbose, ...
                    ['Missing field "', fn, '" added since disgas is not enabled.\n']);
                state.(fn) = 0;
            end
            clear fn
        end
        if model.vapoil
            % RV must be supplied for all cells. This may cause an error.
            model.checkProperty(state, 'rv', nc, 1);
            rvMax = model.getProp(state, 'rvMax');
            [so, sg, rv] = model.getProps(state, 'so', 'sg', 'rv');
            rv(so > 0) = rvMax(so > 0);
            rv = rv.*(sg > 0);
            state = model.setProp(state, 'rv', rv);
        else
            % RS does not really matter. Assign single value.
            fn = model.getVariableField('rv');
            if ~isfield(state, fn)
                dispif(model.verbose, ...
                    ['Missing field "', fn, '" added since vapoil is not enabled.\n']);
                state.(fn) = 0;
            end
            clear fn
        end
    end
    % --------------------------------------------------------------------%
    function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
        [model, state] = prepareTimestep@ReservoirModel(model, state, state0, dt, drivingForces);
    end

    % --------------------------------------------------------------------%
    function [state, report] = updateState(model, state, problem, dx, drivingForces)
        vars = problem.primaryVariables;
        removed = false(size(vars));
        if model.disgas || model.vapoil
            % The VO model is a bit complicated, handle this part
            % explicitly.
            state0 = state;
            state = model.initPropertyContainers(state);

            state = model.updateStateFromIncrement(state, dx, problem, 'pressure', model.dpMaxRel, model.dpMaxAbs);
            state = model.capProperty(state, 'pressure', model.minimumPressure, model.maximumPressure);

            [vars, ix] = model.stripVars(vars, 'pressure');
            removed(~removed) = removed(~removed) | ix;

            % Black oil with dissolution
            [so, sg] = model.getProps(state, 'so', 'sg');
            if model.water
                sw = model.getProp(state, 'sw');
                dsw = model.getIncrement(dx, problem, 'sw');
            else
                sw = 0;
                dsw = 0;
            end
            % Magic status flag, see inside for doc
            st = model.getCellStatusVO(state0, so, sw, sg);

            dr = model.getIncrement(dx, problem, 'x');
            % Interpretation of "gas" phase varies from cell to cell, remove
            % everything that isn't sG updates
            dsg = st{3}.*dr - st{2}.*dsw;

            if model.disgas
                rsMax = model.getProp(state, 'rsMax');
                drs_rel = rsMax.*model.drsMaxRel;
                drs = min(model.drsMaxAbs, drs_rel);
                state = model.updateStateFromIncrement(state, st{1}.*dr, problem, ...
                                                       'rs', inf, drs);
            end

            if model.vapoil
                rvMax = model.getProp(state, 'rvMax');
                drv_rel = rvMax.*model.drsMaxRel;
                drs = min(model.drsMaxAbs, drv_rel);
                state = model.updateStateFromIncrement(state, st{2}.*dr, problem, ...
                                                       'rv', inf, drs);
            end

            dso = -(dsg + dsw);
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

            state = model.updateStateFromIncrement(state, ds, problem, 's', inf, model.dsMaxAbs);
            
            kr = model.FlowPropertyFunctions.RelativePermeability;
            state = kr.applyImmobileChop(model, state, state0);

            % We should *NOT* be solving for oil saturation for this to make sense
            assert(~any(strcmpi(vars, 'so')));
            state = computeFlashBlackOil(state, state0, model, st);
            state.s  = bsxfun(@rdivide, state.s, sum(state.s, 2));

            %  We have explicitly dealt with rs/rv properties, remove from list
            %  meant for autoupdate.
            [vars, ix] = model.stripVars(vars, {'sw', 'so', 'sg', 'rs', 'rv', 'x'});
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
    function scaling = getScalingFactorsCPR(model, problem, names, solver)
        % Get approximate, impes-like pressure scaling factors
        nNames = numel(names);
        
        scaling = cell(nNames, 1);
        handled = false(nNames, 1);
        
        % Take averaged pressure for scaling factors
        state = problem.state;
        fluid = model.fluid;
        if (isprop(solver, 'trueIMPES') || isfield(solver, 'trueIMPES')) && solver.trueIMPES
            % Rigorous pressure equation (requires lots of evaluations)
            rs = state.rs;
            rv = state.rv;
            cfac = 1./(1 - model.disgas*model.vapoil*rs.*rv);
            fp = model.FlowPropertyFunctions;
            b = fp.getProperty(model, state, 'ShrinkageFactors');
            ph = model.getPhaseNames();
            
            iso = ph == 'O';
            isg = ph == 'G';
            isw = ph == 'W';
            [bO, bW, bG] = deal(1);
            if any(isw)
                bW = b{isw};
            end
            if any(iso)
                bO = b{iso};
            end
            if any(isg)
                bG = b{isg};
            end
            sat = problem.state.s;
            rs = rs.*(sat(:, 3) == 0);
            rv = rv.*(sat(:, 2) == 0);
            for iter = 1:nNames
                name = lower(names{iter});
                switch name
                    case 'oil'
                        s = cfac.*(1./bO - model.disgas*rs./bG);
                    case 'water'
                        s = 1./bW;
                    case 'gas'
                        s = cfac.*(1./bG - model.vapoil*rv./bO);
                    otherwise
                        continue
                end
                sub = strcmpi(problem.equationNames, name);
                scaling{iter} = s;
                handled(sub) = true;
            end
        else
            % Very simple scaling factors, uniform over grid
            p = mean(state.pressure);
            useReg = iscell(fluid.bO);
            if useReg
                call = @(x, varargin) x{1}(varargin{:});
            else
                call = @(x, varargin) x(varargin{:});
            end
            for iter = 1:nNames
                name = lower(names{iter});
                switch name
                    case 'oil'
                        if model.disgas
                           rs = call(fluid.rsSat, p);
                           bO = call(fluid.bO,p, rs, true);
                        else
                           bO = call(fluid.bO,p);
                        end
                        s = 1./bO;
                    case 'water'
                        bW = call(fluid.bW,p);
                        s = 1./bW;
                    case 'gas'
                        if model.vapoil
                            rv = call(fluid.rvSat, p);
                            bG = call(fluid.bG, p, rv, true);
                        elseif model.gas
                            bG = call(fluid.bG, p);
                        end
                        s = 1./bG;
                    otherwise
                        continue
                end
                sub = strcmpi(problem.equationNames, name);
                scaling{iter} = s;
                handled(sub) = true;
            end
        end
        if ~all(handled)
            % Get rest of scaling factors from parent class
            other = getScalingFactorsCPR@ReservoirModel(model, problem, names(~handled));
            [scaling{~handled}] = other{:};
        end
    end
    
    function st = getCellStatusVO(model, state, sO, sW, sG)
        status = [];
        if isfield(state, 'status')
            status = state.status;
        end
        st = getCellStatusVO(sO, sW, sG, 'status', status, 'vapoil', ...
                                 model.vapoil, 'disgas', model.disgas);
    end
    
    function [sG, rs, rv, rsSat, rvSat] = calculateHydrocarbonsFromStatusBO(model, ...
                                                          status, sO, x, rs, ...
                                                          rv, pressure)
        [sG, rs, rv, rsSat, rvSat] = calculateHydrocarbonsFromStatusBO(model.fluid, ...
                                                          status, sO, x, rs, ...
                                                          rv, pressure, model.disgas, model.vapoil);
    end
    
    
    function dismat = getDissolutionMatrix(model, rs, rv)
        actPh = model.getActivePhases();
        nPh = nnz(actPh);
        if ~model.disgas
            rs = [];
        end
        if ~model.vapoil
            rv = [];
        end
        
        dismat = cell(1, nPh);
        [dismat{:}] = deal(cell(1, nPh));
        ix = 1;
        jx = 1;
        for i = 1:3
            if ~actPh(i)
                continue
            end
            for j = 1:3
                if ~actPh(j)
                    continue
                end
                if i == 2 && j == 3
                    dismat{i}{j} = rv;
                elseif i == 3 && j == 2
                    dismat{i}{j} = rs;
                end
                jx = jx + 1;
            end
            ix = ix + 1;
        end
    end
    
    function components = getDissolutionMatrixMax(model, pressure)
        [rsMax, rvMax] = deal([]);
        if model.disgas
            rsMax = model.fluid.rsSat(pressure);
        end
        if model.vapoil
            rvMax = model.fluid.rvSat(pressure);
        end
        components = model.getDissolutionMatrix(rsMax, rvMax);
    end
    
end
end

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
