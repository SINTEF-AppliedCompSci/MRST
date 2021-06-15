classdef PressureModel < WrapperModel
    properties
        incTolPressure = 1e-3;
        useIncTol = true;
        reductionStrategy;
        pressureIncTolType = 'relative';
        pressureTol = inf;
        pressureScaling = [];
        laggedStateFunctions = [];
    end
    
    methods
        function model = PressureModel(parent, varargin)
            if isprop(parent, 'useCNVConvergence')
                parent.useCNVConvergence = false;
            end
            model = model@WrapperModel(parent);
            model = merge_options(model, varargin{:});
            if isempty(model.reductionStrategy)
                if isa(parent, 'ThreePhaseBlackOilModel')
                    s = 'analytic';
                else
                    s = 'numerical';
                end
                model.reductionStrategy = s;
            end
            model.AutoDiffBackend = parent.AutoDiffBackend;
        end
        
        function [state, names, origin] = getStateAD(model, state, init)
            if nargin < 3
                init = true;
            end
            % Get the AD state for this model
            [vars, names, origin] = model.getPrimaryVariables(state);
            switch model.reductionStrategy
                case {'analytic', 'numerical'}
                    isP = strcmp(names, 'pressure');
                    origP = origin{isP};
                    keep = isP | cellfun(@(x) ~strcmp(x, origP), origin);
                otherwise
                    error('Unknown reduction strategy %s', model.reductionStrategy);
            end
            pvars = vars(keep);
            pnames = names(keep);
            pscale = model.pressureScaling;
            hasScaling = ~isempty(pscale);
            if init
                if hasScaling
                    % We have a scaling factor to the pressure derivatives.
                    % Scale first, initialize, then divide away. The next fix
                    % comes in the update function, where the scaling is
                    % divided away to keep units consistent.
                    isP = strcmpi(pnames, 'pressure');
                    pvars{isP} = pvars{isP}.*pscale;
                    isBHP = strcmpi(pnames, 'bhp');
                    if any(isBHP)
                        bhpscale = mean(pscale); % To scalar
                        pvars{isBHP} = pvars{isBHP}.*bhpscale;
                    end
                end
                [pvars{:}] = model.AutoDiffBackend.initVariablesAD(pvars{:});
                if hasScaling
                    pvars{1} = pvars{1}./pscale;
                    if any(isBHP)
                        pvars{isBHP} = pvars{isBHP}./bhpscale;
                    end
                end
            end
            vars(keep) = pvars;
            state = model.initStateAD(state, vars, names, origin);
            % Not all were AD-initialized
            names = pnames;
            origin = origin(keep);
            assert(strcmpi(names{1}, 'pressure'));
            if isfield(state, 'laggedStateFunctions')
                % Overwrite initialized contained contents with lagged
                % properties
                fn = fieldnames(state.laggedStateFunctions);
                for i = 1:numel(fn)
                    fld = fn{i};
                    tmp = state.laggedStateFunctions.(fld);
                    sfn = fieldnames(tmp);
                    for j = 1:numel(sfn)
                        n = sfn{j};
                        state.(fld).(n) = state.laggedStateFunctions.(fld).(n);
                    end
                end
            end
        end

        function model = validateModel(model, varargin)
            model = validateModel@WrapperModel(model, varargin{:});
        end

        
        function [eqs, names, types, state] = getModelEquations(pmodel, state0, state, dt, drivingForces)
            model = pmodel.parentModel;
            ncomp = model.getNumberOfComponents();
            [eqs, names, types, state] = model.getModelEquations(state0, state, dt, drivingForces);
            ceqs = eqs(1:ncomp);
            
            w = model.getProp(state, 'PressureReductionFactors');
            % Assemble equations and add in sources
            pressure_equation = 0;
            for i = 1:numel(ceqs)
                pressure_equation = pressure_equation + w{i}.*ceqs{i};
            end
            subs = (ncomp+1):numel(eqs);
            eqs = [{pressure_equation}, eqs(subs)];
            names = ['pressure', names(subs)];
            types = ['cell', types(subs)];
            state = pmodel.updateReductionFactorProps(state);
       end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = getEquations@PhysicalModel(model, state0, state, dt, drivingForces, varargin{:});
        end
        
        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
            [model, state] = prepareTimestep@WrapperModel(model, state, state0, dt, drivingForces);
            state = model.assignReductionFactorProps(state, state0, dt);
            % Ensure that saturations are normalized
            state.s = bsxfun(@rdivide, state.s, sum(state.s, 2));
            % Check for lagged functions - from previous iteration of
            % pressure transport loop, or previous timestep, depending on
            % if the model is linear.
            lfd = model.laggedStateFunctions;
            if ~isempty(lfd)
                % Loop over all groups of lagged state functions and store
                % them in a special field so that they can be set in
                % initStateAD and avoid recomputation
                fn = fieldnames(lfd);
                tmp = state;
                for i = 1:numel(fn)
                    f = fn{i};
                    lagged = lfd.(f);
                    nlag = numel(lagged);
                    if nlag
                        sfg = model.parentModel.(f);
                        fld = sfg.getStateFunctionContainerName();
                        if ~isfield(tmp, fld)
                            tmp = sfg.initStateFunctionContainer(tmp);
                        end
                        for j = 1:nlag
                            [~, tmp] = sfg.get(model.parentModel, tmp, lagged{j});
                        end
                        % Remove other intermediate values that we do not
                        % want to lag! Also, remove caching status.
                        v = struct(tmp.(fld));
                        sfn = fieldnames(v);
                        for j = 1:numel(sfn)
                            if ~any(strcmp(sfn{j}, lagged))
                                v = rmfield(v, sfn{j});
                            end
                        end
                        state.laggedStateFunctions.(fld) = v;
                    end
                end
            end
        end
        
        function state = validateState(model, state)
            state = validateState@WrapperModel(model, state);
            state = model.assignReductionFactorProps(state, state, 1);
        end

        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            pvar = problem.primaryVariables;
            isP = strcmpi(pvar, 'pressure');

            if ~isempty(model.pressureScaling)
                isBHP = strcmp(pvar, 'bhp');
                pscale = model.pressureScaling;
                dx{isP} = dx{isP}./pscale;
                if any(isBHP)
                    dx{isBHP} = dx{isBHP}./mean(pscale);
                end
            end
            p0 = state.pressure;
            [state, report] = model.parentModel.updateState(state, problem, dx, drivingForces);
            if model.verbose > 1
                fprintf('* Pressure: Minimum %f bar, maximum %f bar, mean %f bar\n', min(p0)/barsa, max(p0)/barsa, mean(p0)/barsa)
                sg = model.parentModel.getProp(state, 'sg');
                fprintf('* %d two-phase cells\n', sum(sg > 0));
            end
            if any(isP)
                state.pressureChange = dx{isP};
            else
                state.pressureChange = zeros(size(p0));
            end
        end
        
        function  [convergence, values, names] = checkConvergence(model, problem)
            [convergence, values, names] = model.parentModel.checkConvergence(problem);
            ptol = model.pressureTol;
            itol = model.incTolPressure;
            usePTol = isfinite(ptol);
            useITol = isfinite(itol) && model.useIncTol;
            assert(usePTol || useITol, 'Pressure model has no valid way of checking convergence!');
            if usePTol
                % Check the actual value of the pressure equation. Note:
                % Requires that the pressure reduction factors are scaled
                % relative to each other between nonlinear iterations.
                isP = strcmpi(problem.equationNames, 'pressure');
                [mass, w] = model.parentModel.getProps(problem.state       , ...
                                                 'ComponentTotalMass'      , ...
                                                 'PressureReductionFactors');
                mass  = value(reshape(mass, 1, []));
                w     = value(reshape(w, 1, []));
                scale = problem.dt./sum(mass.*w,2);
                values(1) = norm(value(problem.equations{isP}).*scale, inf);
                convergence(1) = values(1) < ptol;
                names{1} = 'pressure';
            else
                % Skip first equation
                names = names(2:end);
                convergence = convergence(2:end);
                values = values(2:end);
            end
            
            if useITol
                % Check the increment tolerance  for pressure
                if problem.iterationNo > 1 && ~isnan(problem.iterationNo)
                    incType = lower(model.pressureIncTolType);
                    if model.parentModel.G.cells.num == 1 && strcmp(incType, 'relative')
                        % Only norm makes sense
                        incType = 'norm';
                    end
                    dp = problem.state.pressureChange;
                    switch incType
                        case 'relative'
                            p0 = value(problem.state.pressure);
                            range = max(p0) - min(p0);
                            if range == 0
                                range = norm(p0, inf);
                            end
                            scale = range;
                        case 'absolute'
                            % Do nothing
                            scale = 1;
                        case 'norm'
                            scale = norm(value(problem.state.pressure), inf);
                        otherwise
                            error('Unknown pressure increment %s', model.pressureIncTolType);
                    end
                    dp = norm(dp./scale, inf);
                else
                    dp = inf;
                end
                convergence = [dp < itol, convergence];
                names = ['Delta P', names];
                values = [dp, values];
            end
        end

        function [state, report] = updateAfterConvergence(model, varargin)
            [state, report] = updateAfterConvergence@WrapperModel(model, varargin{:});
            % Clean up internal fields
            removeFields = {'laggedStateFunctions', ...
                            'statePressure', ...
                            'pressureChange', ...
                            'reductionFactorProps'};
            for i = 1:numel(removeFields)
                f = removeFields{i};
                if isfield(state, f)
                    state = rmfield(state, f);
                end
            end
            if isfield(state, 'sT')
                state.sT = sum(value(state.s), 2);
            end
            state.statePressure = state;
        end
        
        function rhoS = getSurfaceDensities(model)
            rhoS = model.parentModel.getSurfaceDensities();
        end
        
        function checkPressureReduction(model, state0, state, dt, drivingForces)
            m = model.parentModel;
            state = m.reduceState(state);
            
            state = m.getStateAD(state);
            [M, w] = m.getProps(state, 'ComponentTotalMass', 'PressureReductionFactors');
            M0 = m.getProp(state0, 'ComponentTotalMass');
            e = 0;
            for i = 1:numel(M)
                e = e + (M{i} - M0{i}).*w{i}./dt;
            end
            pos = 1;
            pnames = state.primaryVariables;
            for i = 1:numel(e.jac)
                J = e.jac{i};
                if issparse(J)
                    n = 1;
                    diagonals = {diag(J)};
                    names = pnames{pos};
                    pos = pos + 1;
                else
                    n = J.getNumberOfDiagonals();
                    diagonals = cell(1, n);
                    names = cell(1, n);
                    for j = 1:n
                        diagonals{j} = J.getDiagonalByIndex(j);
                        names{j} = pnames{pos};
                        pos = pos + 1;
                    end
                end
                for j = 1:n
                    d = diagonals{j};
                    name = names{j};
                    if strcmpi(name, 'pressure') && ~isempty(model.pressureScaling)
                        d = d./model.pressureScaling;
                    end
                    if ~isempty(d)
                        fprintf('!! %s: \n * Mean:  %g\n * Min:   %g\n * Max:   %g\n *|Min|:  %g\n *|Mean|: %g\n', ...
                                 name, mean(d), min(d), max(d), min(abs(d)), mean(abs(d)));
                    end
                end
            end
        end
        
        function model = setupStateFunctionGroupings(model, varargin)
            pmodel = model.parentModel;
            pvt = pmodel.PVTPropertyFunctions;
            switch model.reductionStrategy
                case 'analytic'
                    assert(isa(pmodel, 'ThreePhaseBlackOilModel'), ...
                        'Analytical pressure reduction factors currently only implemented for black-oil.');
                case 'numerical'
                    rf = pvt.getStateFunction('PressureReductionFactors');
                    if ~isa(rf, 'NumericalPressureReductionFactors')
                        pvt = pvt.setStateFunction('PressureReductionFactors', NumericalPressureReductionFactors(pmodel));
                    end
                case 'preselected'
                    % Do nothing.
                otherwise
                    error('Unknown reduction strategy');
            end
            model.parentModel.PVTPropertyFunctions = pvt;
        end
        
        function state = assignReductionFactorProps(model, state, state0, dt)
            if strcmpi(model.reductionStrategy, 'numerical')
                mass0 = model.parentModel.getProp(state0, 'ComponentTotalMass');
                props = struct('mass0'   , {mass0}, ...
                               'mass'    , []     , ...
                               'dt'      , dt     , ...
                               'pressure', []     , ...
                               'weights' , []     );
                state.reductionFactorProps = props;
            end
        end

        function state = updateReductionFactorProps(model, state)
            if isfield(state, 'reductionFactorProps')
                pressure = model.getProp(state, 'pressure');
                weights  = model.parentModel.getProp(state, 'PressureReductionFactors');
                mass  = model.parentModel.getProp(state, 'ComponentTotalMass');
                state.reductionFactorProps.pressure = value(pressure);
                state.reductionFactorProps.weights = value(weights');
                state.reductionFactorProps.mass = value(mass');
            end
        end
        
        function model = setLaggedTerms(model, type)
            if isstruct(type)
                lt = type;
            else
                switch lower(type)
                    case 'none'
                        lt = struct('FlowPropertyFunctions', {{}}, ...
                                    'PVTPropertyFunctions', {{}}, ...
                                    'FlowDiscretization', {{}});
                    case 'mobilities'
                        lt = struct('FlowPropertyFunctions', {{'Mobility', 'ComponentMobility'}}, ...
                                    'PVTPropertyFunctions', {{}}, ...
                                    'FlowDiscretization', {{}});
                    otherwise
                        error('Bad type %s', type);
                end
            end
            model.laggedStateFunctions = lt;
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
