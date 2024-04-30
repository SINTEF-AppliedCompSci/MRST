classdef CompositeModel < PhysicalModel
% Model for coupling together multiple submodels, potentially with coupling
% terms between them


    properties
        
        submodels     % Model submodels
        CouplingTerms % Coupling terms between submodels
        
        % Experimental properties for various solution strategies
        solutionStrategy = 'normal';
        solutionSequence = {}
        globalCorrection = true
        nonlinearSolvers
        convergeOuter    = true
        
    end
    
    methods
        %-----------------------------------------------------------------%
        function model = CompositeModel(submodels, varargin)
        % Constructor. Required argument is a cell array of submodels.
        % Coupling terms are currently best set after construction using
        % model.setCouplingTerm.
        
            % Optional input arguments
            opt             = struct('names', {{}});
            [opt, varargin] = merge_options(opt, varargin{:});
            
            % Call parent constructor with empty G and no other arguments
            model = model@PhysicalModel([]);
            if isempty(opt.names)
                % Set submodel names if not already given
                opt.names = cellfun(@class, submodels, 'UniformOutput', false);
                % Check that all submodel names are unique
                assert(numel(opt.names) == numel(unique(opt.names)), ...
                    ['Repeated submodel names, please provide unique ', ...
                     'names using optional input `names`']);
            end
            
            % Set submodels to model
            model.submodels = struct();
            for i = 1:numel(submodels)
                model.submodels.(opt.names{i}) = submodels{i};
            end
            
            % Set up state function grouping for coupling terms
            model.CouplingTerms = StateFunctionGrouping('CouplingTerms');
            
            % Set autodiff backend from first submodel
            mnames = model.getModelNames();
            model.AutoDiffBackend = model.submodels.(mnames{1}).AutoDiffBackend;
            
            % Merge options
            model = merge_options(model, varargin{:});

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function model = setCouplingTerm(model, coupling, name)
        % Set coupling terms to the model
        
            if iscell(coupling)
                % Cell array of coupling terms, set each of them
                for i = 1:numel(coupling)
                    model = model.setCouplingTerm(coupling{i}, name{i});
                end
                return;
            end
            
            % Set name of coupling term as classname if not given
            if nargin < 3, name = class(coupling); end
            % Check that name has not already been used
            names = model.CouplingTerms.getNamesOfStateFunctions();
            assert(~ismember(name, names), ...
                ['A coupling term with name %s already exists. ', ...
                 'Please provide a unique name with as a third ', ...
                 'argument to this function'], name);
             
            % Set to coupling term state function grouping
            model.CouplingTerms ...
                = model.CouplingTerms.setStateFunction(name, coupling);
            
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        % Driving forces
        %-----------------------------------------------------------------%
        function forces = validateDrivingForces(model, forces, index)
        
            forces = validateDrivingForces@PhysicalModel(model, forces, index);
            names = model.getModelNames();
            for name = names
                forces.(name{1}) = model.submodels.(name{1}).validateDrivingForces(forces.(name{1}), index);
            end
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function forces = getValidDrivingForces(model)
        % Validate all driving forces
        
            forces = getValidDrivingForces@PhysicalModel(model);
            names = model.getModelNames();
            for name = names
                forces.(name{1}) = model.submodels.(name{1}).getValidDrivingForces();
            end
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [model, state] = updateForChangedControls(model, state, forces)
        % Update for changed controls during a simulation
            
            names = model.getModelNames();
            for name = names
                [model.submodels.(name{1}), state.(name{1})] ...
                    = model.submodels.(name{1}). ...
                    updateForChangedControls(state.(name{1}), forces.(name{1}));
            end
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        % State
        %-----------------------------------------------------------------%
        function state = validateState(model, state)
        % Validate state and check if it is ready for simulation
            
            names = model.getModelNames();
            for name = names
                state.(name{1}) = model.submodels.(name{1}).validateState(state.(name{1}));
            end
            state = validateState@PhysicalModel(model, state);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [state, names, origin] = getStateAD(model, state, init)
        % Get the AD state for all submodels, with cross-model derivatives
        
            if nargin < 3, init = true; end
            
            active = model.getActiveModels(state);
            nm     = numel(active);
            [vars, names, origin] = deal(cell(1,nm));
            
            for i = 1:nm
                % Get the AD state for this model
                mname = active{i};
                [vars{i}, names{i}, origin{i}] ...
                    = model.submodels.(mname).getPrimaryVariables(state.(mname));
            end

            if init
                % Assemble all primary variables into one vector and
                % initialize as AD
                [varsAD, mapl, mapn] = unpackVars(vars, 1, 1, [], []);
                [varsAD{:}] = model.AutoDiffBackend.initVariablesAD(varsAD{:});
                vars = packVars(varsAD, mapl, mapn, 2);
                state = model.setPrimaryVariables(state, names);
            end
            
            % Initialize AD state for all submodels
            state = model.initStateAD(state, vars, names, origin);
            
            % Evaluate all statefunctions for all models if requested
            if strcmpi(model.stateFunctionEvaluationMode, 'full')
                state = model.evaluateAllStateFunctions(state);
            end

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [vars, names, origin] = getPrimaryVariables(model, state)
            
            active = model.getActiveModels(state);
            nm     = numel(active);
            [vars, names, origin] = deal(cell(1,nm));
            for i = 1:nm
                % Get the AD state for this model
                mname = active{i};
                [vars{i}, names{i}, origin{i}] ...
                    = model.submodels.(mname).getPrimaryVariables(state.(mname));
            end
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function state = setPrimaryVariables(model, state, names)
            
            active = model.getActiveModels(state);
            nm     = numel(active);
            for i = 1:nm
                mname = active{i};
                if isa(model.submodels.(mname), 'CompositeModel')
                    state.(mname) = model.submodels.(mname).setPrimaryVariables( ...
                        state.(mname), names{i});
                else
                    state.(mname).primaryVariables = names{i};
                end
            end
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function state = initStateAD(model, state, vars, names, origin)
            
            active = model.getActiveModels(state);
            nm     = numel(active);
            for i = 1:nm
                mname = active{i};
                state.(mname) = model.submodels.(mname).initStateAD( ...
                    state.(mname), vars{i}, names{i}, origin{i});
            end
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function state = evaluateAllStateFunctions(model, state)
            
            active = model.getActiveModels(state);
            nm = numel(active);
            for i = 1:nm
                mname = active{i};
                state.(mname) = model.submodels.(mname). ...
                        evaluateAllStateFunctions(state.(mname));
            end
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%       
        function state = reduceState(model, state, removeContainers)
        % Reduce state to doubles, and optionally remove the property
        % containers to reduce storage space
        
            % Reduce each submodel state
            for mname = model.getModelNames()
                state.(mname{1}) = model.submodels.(mname{1}).reduceState( ...
                    state.(mname{1}), removeContainers);
            end
            % Reduce state itself
            state = reduceState@PhysicalModel(model, state, removeContainers);
            
        end    
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, forces)
        % Update state based on linearized update
            
            [problems, map] = model.splitProblem(problem);
            dxs = cell(1, model.numModels);
            active = model.getActiveModels(problem.state);
            nm = numel(active);
            for i = 1:nm
                ix = map == i;
                dxs{i} = dx(ix);
            end
            
            report = struct();
            for i = 1:nm
                mname = active{i};
                [state.(mname), report.(mname)] ...
                    = model.submodels.(mname).updateState( ...
                    state.(mname), problems{i}, dxs{i}, forces.(mname));
            end
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [state, report] = stepFunction(model, state, state0, dt, drivingForces, linsolver, nonlinsolver, iteration, varargin)
        % Perform one nonlinear update. If model.solutionStrategy is set to
        % `normal`, the function performs a standard Newton update, with
        % damping as prescribed by each submodel.
        %
        % Experimental features: If set to `sequential`, the model will
        % solve the entire timestep for each submodel sequentially, in the
        % order prescribed by model.solutionSequence. If set to
        % `field-split`, the model will do a full global update after
        % solving the entire timestep for each submodel.
        
            if strcmpi(model.solutionStrategy, 'normal') || isfield(state, 'activeModels')
                [state, report] = stepFunction@PhysicalModel( ...
                    model, state, state0, dt, drivingForces, ...
                    linsolver, nonlinsolver, iteration, varargin{:} ...
                );
                return
            end
            
            % Experimental solution strategies (sequential and field-split)
            is_seq = strcmpi(model.solutionStrategy, 'sequential');
            is_fs  = strcmpi(model.solutionStrategy, 'field-split');
            assert(any([is_seq, is_fs]));
            
            outOfIterations = iteration > nonlinsolver.maxIterations;
            % Check convergence
            problem = model.getEquations(state0, state, dt, drivingForces, 'resOnly', true, varargin{:});
            [convergence, values, resnames] = model.checkConvergence(problem);

            report = model.makeStepReport(...
                        'Residuals',    values, ...
                        'ResidualsConverged', convergence);
            
            % Minimum number of iterations can be prescribed, i.e., we
            % always want at least one set of updates regardless of
            % convergence criterion.
            doneMinIts = iteration > nonlinsolver.minIterations;
            
            doGlobal = model.globalCorrection || is_fs;
            
            converged = false(numel(model.solutionSequence),1);
            
            if (~(all(convergence) && doneMinIts) && ~outOfIterations)

                doLocal = is_seq || (is_fs && iteration == 1);
                if doLocal
                    for i = 1:numel(model.solutionSequence)

                        mname = model.solutionSequence{i};
                        state.activeModels = mname;

                        forceArg = getDrivingForces(model, drivingForces);
                        
                        sname = strjoin(mname, '_');
                        nls = model.nonlinearSolvers.(sname);

                        model.verbose = false;
                        [state, report.(mname{1})] = nls.solveTimestep(state0, dt, model, 'initialGuess', state, forceArg{:});
                        model.verbose = true;
                        state = model.reduceState(state, true);
                        state = rmfield(state, 'activeModels');

                        converged(i) = report.(mname{1}).Converged;
                        if ~converged(i)
                            break;
                        else
%                             model.stepFunctionIsLinear = true;
                        end

                    end
                    if is_seq, doGlobal = doGlobal && all(converged); end
                end
                
                if doGlobal
                    state.Reservoir = model.submodels.Reservoir.initStateFunctionContainers(state.Reservoir);
                    state.Wellbore  = model.submodels.Wellbore.initStateFunctionContainers(state.Wellbore);
                    [state, report] = stepFunction@PhysicalModel(model, state, state0, dt, drivingForces, linsolver, nonlinsolver, iteration, varargin{:});
                else
                    report.Failure = ~all(converged);
                end
                state = model.reduceState(state, true);
                
            end

            modelConverged = all(convergence);
            isConverged = (modelConverged && doneMinIts) || model.stepFunctionIsLinear ...
                || (all(converged) && ~model.convergeOuter);
            
            report.Converged = isConverged;
            
            if model.verbose && (~doGlobal || isConverged)
                printConvergenceReport(resnames, values, convergence, iteration, isConverged);
            end
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        % Equations
        %-----------------------------------------------------------------%        
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces)
        % Get model equations for each submodel, with any prescribed
        % coupling terms inserted
            
            % Get equations for each active model
            active = model.getActiveModels(state);
            nm     = numel(active);
            [eqs, names, types] = deal(cell(1, nm));
            for i = 1:nm
                mname = active{i};
                [eqs{i}, names{i}, types{i}, state.(mname)] ...
                    = model.submodels.(mname).getModelEquations( ...
                    state0.(mname), state.(mname), dt, drivingForces.(mname));
            end
            
            % Insert any couplings between models
            cnames = model.CouplingTerms.getNamesOfStateFunctions();
            if isempty(cnames), return; end
            for cname = cnames'
                % Evaluate coupling term state function
                [q, state] = model.getProp(state, cname{1});
                % Insert into correct equation
                ct  = model.CouplingTerms.getStateFunction(cname{1});
                eqs = ct.insertCoupling(model, eqs, active, names, q);
            end
            
        end
        %-----------------------------------------------------------------%
                
        %-----------------------------------------------------------------%
        function problem = setupLinearizedProblem(model, eqs, types, names, primaryVars, state, dt)
        % Set up linearized problem with all equations and primary
        % variables of all submodels
        
            % Get names of active models and initialize
            active   = model.getActiveModels(state);
            nm       = numel(active);
            problems = cell(nm, 1);
            
            % Set up linearized problem for each submodel
            for i = 1:nm
                mname = active{i};
                problems{i} = model.submodels.(mname).setupLinearizedProblem( ...
                    eqs{i}, types{i}, names{i}, primaryVars{i}, state.(mname), dt);
            end
            
            % Gather into composite linearized problem
            problem  = LinearizedProblem({}, []);
            pressure = [];
            for i = 1:nm
                eqNames = problems{i}.equationNames;
                eqNames = cellfun(@(names) [active{i}, '.', names], eqNames, 'UniformOutput', false);                
                problem = problem.appendEquations(problems{i}.equations, problems{i}.types, eqNames);
                problem.primaryVariables = [problem.primaryVariables, problems{i}.primaryVariables];
                if isfield(problems{i}.state, 'pressure')
                    % Stack pressures of each submodel for use in CPR
                    pressure = [pressure; problems{i}.state.pressure]; %#ok
                end
                problems{i} = [];
            end
            % Set problem properties
            problem.dt             = dt;
            problem.state          = state;
            problem.state.pressure = pressure; 
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [convergence, values, names] = checkConvergence(model, problem, varargin)
        % Check convergence of composite problem
            
            % Split into subproblems
            problems = model.splitProblem(problem);
            
            % Allocate
            mnames = model.getActiveModels(problem.state);
            nm = numel(mnames);
            [convergence, values, names] = deal(cell(1,nm));
            for i = 1:nm
                mname = mnames{i};
                [convergence{i}, values{i}, names{i}] ...
                    = model.submodels.(mname).checkConvergence(problems{i}, varargin{:});
                names{i} = cellfun(@(name) [mname, '.', name], names{i}, 'UniformOutput', false);
            end
            convergence = horzcat(convergence{:});
            values      = horzcat(values{:});
            names       = horzcat(names{:});

        end    
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [problems, map] = splitProblem(model, problem, removeCouplings)
        % Split composite problem into subproblems
            
            if nargin < 3, removeCouplings = false; end
            mnames = model.getActiveModels(problem.state);
            nm = numel(mnames);
            problems = cell(1, nm);
            map = nan(problem.numeq, 1);
            for i = 1:nm
                mname = mnames{i};
                startswith = @(str, x) strcmp(str(1:min(numel(str), numel(x)+1)), [x, '.']);
                ix = find(cellfun(@(x) startswith(x, mname), problem.equationNames));
%                 names = strrep(problem.equationNames(ix), [mname, '.'], '');
                names = cellfun(@(name) name(numel([mname, '.'])+1:end), problem.equationNames(ix), 'UniformOutput', false);
                problems{i} = LinearizedProblem( ...
                    problem.equations(ix)       , ...
                    problem.types(ix)           , ...
                    names                       , ...
                    problem.primaryVariables(ix), ...
                    problem.state.(mname)       , ...
                    problem.dt                    ...
                );
                map(ix) = i;
                if removeCouplings
                    eqs = problems{i}.equations;
                    for j = 1:numel(eqs)
                        eqs{j}.jac = eqs{j}.jac(ix);
                    end
                    problems{i}.equations = eqs;
                end
            end
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [p, state] = getProp(model, state, name)
        % Get a single property from the nonlinear state
            
            name = regexp(name, '[.]', 'split');
            if numel(name) >= 2
                % We are asking for a submodel property, query the submodel
                [p, state.(name{2})] = model.submodels.(name{1}).getProp(state.(name{1}), name{2});
            else
                % Property stems from this model (likely a coupling term)
                [p, state] = getProp@PhysicalModel(model, state, name{1});
            end
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [groups, names, models] = getStateFunctionGroupings(model)
        % Return model state function groupings
        
            groups = {model.CouplingTerms};
            names  = {'CouplingTerms'};
            models = {model};
            
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function model = validateModel(model, varargin)
        % Validate model and all submodels
        
            hasForces = nargin > 1;
            if hasForces
                forces = varargin{1};
                hasForces = true;
                varargin = varargin(2:end);
            end
            names = model.getModelNames();
            for name = names
                args = varargin;
                if hasForces, args = [{forces.(name{1})}, varargin]; end
                model.submodels.(name{1}) = model.submodels.(name{1}).validateModel(args{:});
            end
            
        end
        %-----------------------------------------------------------------%

        % ----------------------------------------------------------------%
        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
        % Prepare timestep for all submodels
            
            mnames = model.getModelNames();
            for mname = mnames
                [model.submodels.(mname{1}), state.(mname{1})] ...
                    = model.submodels.(mname{1}).prepareTimestep( ...
                    state.(mname{1}), state0.(mname{1}), dt, drivingForces.(mname{1}));
            end
            
        end
        % ----------------------------------------------------------------%
        
        
        %-----------------------------------------------------------------%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
         % Final update to the state after convergence has been achieved
         
            % Update according to each submodel
            mnames = model.getModelNames();
            for mname = mnames
                [state.(mname{1}), report.(mname{1})] ...
                    = model.submodels.(mname{1}).updateAfterConvergence( ...
                    state0.(mname{1}), state.(mname{1}), dt, drivingForces.(mname{1}));
            end
            % Multimodel update
            state = updateAfterConvergence@PhysicalModel(model, state0, state, dt, drivingForces);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [lambda, lambdaVec, report] = solveAdjoint(model, solver, getState, getObj, schedule, gradient, stepNo, varargin)
            
            validforces = model.getValidDrivingForces();
            dt_steps = schedule.step.val;
            
            current = getState(stepNo);
            before  = getState(stepNo - 1);
            dt = dt_steps(stepNo);
            
            lookupCtrl = @(step) schedule.control(schedule.step.control(step));
            % get forces and merge with valid forces
            forces = model.getDrivingForces(lookupCtrl(stepNo));
            forces = merge_options(validforces, forces{:});
            
            % search for reservoir-models
            modNames    = fieldnames(model.submodels);
            hasFacility = cellfun(@(nm)isa(model.submodels.(nm).internalModel, 'ReservoarModel') ...
                                    && ~isempty(model.submodels.(nm).internalModel.FacilityModel), modNames);
            for km = 1:numel(modNames)
                frc = forces.(modNames{km});
                if hasFacility(km)
                    model.submodels.(modNames{km}).internalModel.FacilityModel = ...
                        model.submodels.(modNames{km}).internalModel.FacilityModel.validateModel(frc);
                else
                    model.submodels.(modNames{km}) = model.submodels.(modNames{km}).validateModel(frc);
                end
            end
            
            % Initial state typically lacks wellSol-field, so add if needed
            if stepNo == 1
                before = model.validateState(before);
            end
            
            problem = model.getAdjointEquations(before, current, dt, forces, ...
                'reverseMode', false, 'iteration', inf);
            
            if stepNo < numel(dt_steps)
                after   = getState(stepNo + 1);
                dt_next = dt_steps(stepNo + 1);
                current = model.getReverseStateAD(current);
                %assert(isequal(current.primaryVariables, primaryVars))
                % check if controls for for stepNo/stepNo+1 are different
                if diff(schedule.step.control([stepNo;stepNo+1])) == 0
                    forces_p = forces;
                else
                    % get forces and merge with valid forces
                    forces_p = model.getDrivingForces(lookupCtrl(stepNo + 1));
                    forces_p = merge_options(validforces, forces_p{:});
                    
                    for km = 1:numel(modNames)
                        frc_p = forces_p.(modNames{km});
                        if hasFacility(km)
                            model.submodels.(modNames{km}).internalModel.FacilityModel = ...
                                model.submodels.(modNames{km}).internalModel.FacilityModel.validateModel(frc_p);
                        else
                            model.submodels.(modNames{km}) = model.submodels.(modNames{km}).validateModel(frc_p);
                        end
                    end
                end
                problem_p = model.getAdjointEquations(current, after, dt_next, forces_p,...
                    'iteration', inf, 'reverseMode', true);
            else
                problem_p = [];
            end
            [lambda, lambdaVec, rep] = solver.solveAdjointProblem(problem_p,...
                problem, lambda, getObj(stepNo,model,problem.state), model, ...
                varargin{:});
            report = struct();
            report.Types = problem.types;
            report.LinearSolverReport = rep;
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [problem, state] = getAdjointEquations(model, state0, state, dt, forces, varargin)
            
            opt = struct('Verbose',     mrstVerbose,...
                         'reverseMode', false,...
                         'resOnly',     false, ...
                         'iteration',   -1);
            opt = merge_options(opt, varargin{:});
            % assume primary vars do not change
            if opt.reverseMode
                [state, primaryVars]  = model.getStateAD(state, false);
            else
                [state, primaryVars] = model.getStateAD(state, ~opt.resOnly);
            end
            [eqs, names, types, state] = model.getModelEquations(state0, state, dt, forces);
            problem = model.setupLinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
        end     
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        % Useful utility functions
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function n = numModels(model)
        % Get number of submodels
        
            names = model.getModelNames();
            n = numel(names);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function names = getModelNames(model)
        % Get cell array of submodel names
            
            names = fieldnames(model.submodels)';
            
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function names = getActiveModels(model, state)
        % Get active submodels based on `activeModels` field on state. If
        % no such field is present, all models are assumed to be active
            
            names = model.getModelNames();
            if isfield(state, 'activeModels')
                active = ismember(names, state.activeModels);
                names = names(active);
            end
            
        end
        %-----------------------------------------------------------------%
        
    end
    
end

%-------------------------------------------------------------------------%
function [vars, mapl, mapn] = unpackVars(vars, level, num, mapl, mapn)

    vars0 = vars;
    vars = {};
    for i = 1:numel(vars0)
        if iscell(vars0{i})
            [v, mapl, mapn] = unpackVars(vars0{i}, level+1, i, mapl, mapn);
            vars = [vars, v];
        else
            vars = [vars, vars0{i}];
            mapl = [mapl, level];
            mapn = [mapn, num];
        end
    end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function vars = packVars(varsAD, mapl, mapn, level)

    subsetl = mapl == level;
    n = max(mapn(subsetl));
    vars = cell(1,n);
    vl = varsAD(subsetl);
    for node = 1:n
        subsetn = mapn(subsetl) == node;
        if any(subsetn)
            vn = vl(subsetn);
        else
            vn = packVars(varsAD(~subsetl), mapl(~subsetl), mapn(~subsetl), level+1);
        end
        vars{node} = vn;
    end

end
%-------------------------------------------------------------------------%

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
