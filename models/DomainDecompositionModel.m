classdef DomainDecompositionModel < WrapperModel
    % Nonlinear domain decomposition (NLDD) solver for generic MRST models
    properties
        strategy = 'additive' % Additive or multiplicative NLDD
        parallel = false      % Solve subdomains concurrently (additive only)
        partition             % Partition defining subdomains (instance of class Partition)
        subdomainSetup        % Cell array of subdomain setups
        submodelFn            % Function handle to set submodel. Empty will default to using SubdomainModel
        % Subdomain properties
        overlap         = 0  % Overlap in cells into neighboring subdomains
        verboseSubmodel = 0  % Verbosity of subdomain model solves:
                             % 0 <-> none, 1 <-> nonlinear solver, 2 <-> model & nonlinear solver
        subdomainTol    = [] % Subdomain tolerances relative to parent
                             % model. Given as cell array on the form
                             % {'name', fraction}, e.g., {'toleranceCNV', 0.1};
        useGlobalPressureRange = false % Use pressure range of entire mode
                                       % for models with relative presure
                                       % tolerances
        % Parallel properties
        pool       % Paralell pool
        domains    % Domains per worker
        localModel % Local copy of model (one for each worker)
        NumWorkers % Number of workers
    end
    
    methods
        %-----------------------------------------------------------------%
        function model = DomainDecompositionModel(parent, partition, varargin)
            % Base model initialization
            model = model@WrapperModel(parent);
            % Set partition
            if isnumeric(partition)
                % Construct static Partition from partition vector
                assert(numel(partition) == parent.G.cells.num, ...
                    'Numerical input argument ''partition'' must have G.cells.num entries.')
                partition = StaticPartition(partition);
            else
                % Check class
                assert(isa(partition, 'Partition'));
            end
            model.partition = partition;
            % Merge options
            model = merge_options(model, varargin{:});
            % Control input
            assert(any(strcmpi(model.strategy, {'additive', 'multiplicative'})), ...
                'Supported strategies are ''additive'' and ''multiplicative''');
            if model.parallel
                % Check if parallel toolbox is avialiable
                ok = ~isempty(ver('parallel'));
                if ~ok, msg = 'Parallel computing toolbox not aviailable.'; end
                % Check strategy
                ok = ok & ~strcmpi(model.strategy, 'multiplicative');
                if ~ok, msg = 'Parallel implemetation does not support multiplicative strategy.'; end
                % Set parallel
                if ~ok, warning([msg, ' Switching to serial']); end
                model.parallel = ok;
            end
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = stepFunction(model, state, state0, dt, drivingForces, linsolver, nonlinsolver, iteration, varargin)
            % Prepare state
            state      = model.prepareState(state, iteration);
            iterations = state.iterations;
            % Get residual
            outOfIterations = iteration > nonlinsolver.maxIterations;
            % Check convergence
            [convergence, values, resnames, state] = model.convergenceCheck(state0, state, dt, drivingForces, iteration, varargin{:});
            % Minimum number of iterations can be prescribed, i.e., we
            % always want at least one set of updates regardless of
            % convergence criterion.
            doneMinIts = iteration > nonlinsolver.minIterations;
            subdomainReport = [];
            if (~(all(convergence) && doneMinIts) && ~outOfIterations)
                % Solve subdomains
                timer = tic();
                [state, subdomainReport] = model.solveSubDomains(state0, dt, drivingForces, state);
                t_solve = toc(timer);
                subdomainReport.WallTime = t_solve;
                state.iterations = iterations + subdomainReport.Iterations;
            end
            modelConverged = all(convergence);
            isConverged = (modelConverged && doneMinIts) || model.stepFunctionIsLinear;
            % If step function is linear, we need to call a residual-only
            % equation assembly to ensure that indirect/derived quantities are
            % set with the updated values (fluxes, mobilities and so on).
            if model.stepFunctionIsLinear
                [~, state] = model.getEquations(state0, state, dt, drivingForces, ...
                                       'ResOnly', true, ...
                                       'iteration', iteration+1, ...
                                       varargin{:});
               state = model.reduceState(state, true);
            end
            if model.verbose
                printConvergenceReport(resnames, values, convergence, iteration, isConverged);
            end
            % Make report
            report = model.makeStepReport('Converged'           , isConverged , ...
                                          'Solved'              , ~isConverged, ...
                                          'ResidualsConverged'  , convergence , ...
                                          'Residuals'           , values      );
            report.SubdomainReport = subdomainReport;
        end
        
        %-----------------------------------------------------------------%
        function [convergence, values, resnames, state, problem] = convergenceCheck(model, state0, state, dt, drivingForces, iteration, varargin)
            % Check convergence of the full problem. Note that we evaluate
            % the residuals without Jacobians
            problem = model.getEquations(state0, state, dt, drivingForces, ...
                                         'ResOnly'   , true              , ...
                                         'iteration' , iteration         , ...
                                          varargin{:}                    );
            problem.iterationNo   = iteration;
            problem.drivingForces = drivingForces;
            % Check convergence
            [convergence, values, resnames] = model.parentModel.checkConvergence(problem);
        end
        
        %-----------------------------------------------------------------%
        function state = prepareState(model, state, iteration)
            % Initialize total subdomain iterations
            if iteration == 1
                state.iterations = zeros(model.G.cells.num,1);
            end
            % Remove fractional derivatives if they exist
            if isfield(state, 'FractionalDerivatives')
                state = rmfield(state, 'FractionalDerivatives');
            end
            % Remove statePressure
            if isfield(state, 'statePressure')
                state = rmfield(state, 'statePressure');
            end
            
            % Use global pressure range when computing dpRel
            if model.useGlobalPressureRange
                range = max(state.pressure) - min(state.pressure);
                rmodel = getReservoirModel(model);
                tol = inf;
                if isprop(rmodel, 'incTolPressure')
                    tol   = rmodel.incTolPressure;
                end
                if isinf(tol)
                    tol = 1e-3;
                end
                range = max(range, mean(state.pressure)*tol);
                state.pressureRange = range;
            end
%             state = stripState(state);
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = solveSubDomains(model, state0, dt, drivingForces, state)
            % Solve the subdomains
            if ~model.parallel
                [state, report] = model.solveSubDomainsSerial(state0, dt, drivingForces, state);
            else
                [state, report] = model.solveSubDomainsParallel(state0, dt, drivingForces, state);
            end
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = solveSubDomainsSerial(model, state0, dt, drivingForces, state)
            % Solve the subdomains in serial mode
            % Initialize
            iterations = nan(model.G.cells.num,1);
            [stateInit, stateFinal] = deal(state);
            subreports = cell(max(model.partition.value), 1);
            for i = 1:max(model.partition.value)
                % Solve subdomain
                [stateInit, stateFinal, subreports{i}, iterations] = model.solveSubDomain(model.subdomainSetup{i}, state0, dt, drivingForces, stateInit, stateFinal, iterations);
            end
            state = stateFinal;
            % Make step report
            report = model.makeSubdomainStepReport('Iterations', iterations  , ...
                                                   'subreports', {subreports});
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = solveSubDomainsParallel(model, state0, dt, drivingForces, state)
            % Solve the subdomains in parallel mode using spmd
            % Initialize
            stateInit = stripState(state);
            state0    = stripState(state0);
            sds = model.subdomainSetup;
            lm  = model.localModel;
            add = model.getAddStatesFun(stateInit);
            ticBytes(model.pool);
            spmd
                localTimer = tic(); % Time the local solve (db only)
                map        = [];
                stateFinal = stateInit;
                iterations = zeros(lm.G.cells.num,1);
                subreports = cell(numel(lm.domains),1);
                for i = 1:numel(lm.domains)
                    % Solve subdomain
                    [stateInit, stateFinal, subreports{i}, iterations, submodel] = lm.solveSubDomain(sds{i}, state0, dt, drivingForces, stateInit, stateFinal, iterations);
                    map = mergeMappings(map, submodel.mappings);
                end
                t_localSolve = toc(localTimer);
                % Sum iterations and gather on worker 1
                iterations = gop(@plus, iterations, 1);
                if numel(lm.domains) > 0
                    stateFinal = getSubState(stateFinal, map);
                end
                in = struct('mappings', map, 'state', stateFinal);
                out = gop(add, in, 1);
                if labindex == 1
                    stateFinal = out.state;
                end
                t_localComm = toc(localTimer) - t_localSolve; %#ok
            end
            bytes = tocBytes(model.pool); %#ok
            timer  = tic();
            iterations = iterations{1};
            state      = stateFinal{1};
            t_comm     = toc(timer); %#ok
            % Make step report
            subreports = subreports(:);
            subreports = vertcat(subreports{:});
            report = model.makeSubdomainStepReport('Iterations', iterations  , ...
                                                   'subreports', {subreports});
        end
        
        %-----------------------------------------------------------------%
        function [stateInit, stateFinal, subreport, iterations, varargout] = solveSubDomain(model, setup, state0, dt, drivingForces, stateInit, stateFinal, iterations)
            % Solve a single subdomain
            % Get submodel
            if isempty(setup.Model)
                setup = model.getSubdomainSetup(setup.Number, true);
            end
            submodel = setup.Model;
            mappings = submodel.mappings;
            % Get subdomain forces
            [subforces, mappings] = getSubForces(drivingForces, mappings);
            % Set facility model
            submodel = model.setSubFacilityModel(submodel, mappings);
            % Get subdomain states
            [substateInit, mappings] = getSubState(stateInit, mappings);
            substate0                = getSubState(state0   , mappings);
            % Solve timestep
            nls      = setup.NonlinearSolver;
            forceArg = getDrivingForces(submodel, subforces);
            [substate, subreport] = nls.solveTimestep(substate0, dt, submodel, 'initialGuess', substateInit, forceArg{:});
            if ~subreport.Converged
               substate = substateInit; 
            end
            % Map substate to state
            state = mapState(stateInit, substate, mappings);
            % Update initial and final state
            switch model.strategy
                case 'additive'
                    % Additive method: map substate to final state
                    if subreport.Converged
                        stateFinal = mapState(stateFinal, substate, mappings);
                    end
                case 'multiplicative'
                    % Multiplicative: update initial and final state
                    [stateInit, stateFinal] = deal(state);
            end
            % Update iterations
            cells = model.partition.value == setup.Number;
            iterations(cells) = subreport.Iterations;
            % Store number of subdomain cells in subreport
            subreport.nc = nnz(mappings.cells.internal);
            % Extra output if requested
            varargout = cell(1,nargout-3);
            if nargout > 3
                submodel.mappings = mappings;
                varargout{1} = submodel;
                varargout{2} = state;
                varargout{3} = substate0;
                varargout{4} = subforces;
            end
        end
        
        %-----------------------------------------------------------------%
        function report = makeSubdomainStepReport(model, varargin)%#ok
            % Make subdomain step report
            report = struct('Iterations', [], 'subreports', [], 'WallTime', []);
            report = merge_options(report, varargin{:});
        end
        
        %-----------------------------------------------------------------%
        function model = updateSubdomainSetupSerial(model)
            % Update subdomain setups in serial mode
            p = model.partition.value;
            setup = cell(max(p),1);
            for i = 1:max(p)
                setup{i} = model.getSubdomainSetup(i);
            end
            model.subdomainSetup = setup;
        end
        
        %-----------------------------------------------------------------%
        function model = updateSubdomainSetupParallel(model)
            % Update subdomain setups in parallel mode
            p  = model.partition;
            lm = model.localModel;
            spmd
                lm.partition = p;
                dom = getLocalDomains(lm.NumWorkers, lm.partition.value, labindex);
                setup = cell(1,numel(dom));
                for i = 1:numel(dom)
                    setup{i} = lm.getSubdomainSetup(dom(i));
                    G = setup{i}.Model.G;
                    G = rmfield(G, 'nodes');
                    setup{i}.Model.G = G;
                    setup{i}.Model.parentModel.G = G;
                end
                lm.domains = dom;
            end
            model.localModel = lm;
            model.subdomainSetup = setup;
        end
        
        %-----------------------------------------------------------------%
        function setup = getSubdomainSetup(model, i, compute)
            % Get subdomain setup with simulation-ready submodel and
            % suitable linear and nonlinear solver
            setup = struct('Model', [], 'NonlinearSolver', [], 'Number', i);
            p = model.partition.value;
            if nargin < 3
                compute = false;
            end
            if max(p) > 500 && ~compute
                % Don't precomute setups if there are more than 500 of them
                return
            end
            % Make submodel
            verbose  = model.verboseSubmodel; % Subdomain solve verbosity
            cells    = p == i;                % Cell subset
            if isempty(model.submodelFn)
                submodel = SubdomainModel(model.parentModel, cells, ...
                                          'overlap', model.overlap, ...
                                          'verbose', verbose == 2);
            else
                submodel = model.submodelFn(model, cells);
            end
            if isa(submodel.parentModel, 'SubdomainModel')
                % Update mapping externals if we parent is also a submodel
                pmappings = submodel.parentModel.mappings;
                % Cell mappings
                cex = pmappings.external(pmappings.cells.keep) & submodel.mappings.cells.keep;
                submodel.mappings.cells.external = submodel.mappings.cells.external | cex;
                % Face mappings
                fex = pmappings.external(pmappings.faces.keep) & submodel.mappings.faces.keep;
                submodel.mappings.faces.external = submodel.mappings.faces.external | fex;
            end
            % Adjust submodel tolerances
            submodel = model.setSubdomainTolerances(submodel);
            % Get linear solver
            rmodel = getReservoirModel(submodel);
            ncomp  = rmodel.getNumberOfComponents();
            ndof   = rmodel.G.cells.num*ncomp;
            if isa(submodel.parentModel, 'ReservoirModel')
                lsol = selectLinearSolverAD(rmodel, 'verbose', verbose == 1);
            elseif isa(submodel.parentModel, 'PressureModel')
                if ndof > 1e4
                    lsol = AMGCLSolverAD('tolerance', 1e-4, 'verbose', verbose == 1);
                else
                    lsol = BackslashSolverAD();
                end
            else
                lsol = BackslashSolverAD();
            end
            isTransport = isa(submodel.parentModel, 'TransportModel');
            % Get nonlinear solver
            nls = NonLinearSolver('minIterations' , 0                         , ...
                                  'maxIterations' , 25                        , ...
                                  'ErrorOnFailure', false                     , ...
                                  'verbose'       , verbose                   , ...
                                  'useRelaxation' , isTransport               , ...
                                  'LinearSolver'  , lsol                      , ...
                                  'identifier'    , ['SUBDOMAIN ', num2str(i)]);
            % Make subdomain setup
            setup.Model = submodel;
            setup.NonlinearSolver = nls;
        end
        
        %-----------------------------------------------------------------%
        function model = validateModel(model, varargin)
            % Validate model
            model = validateModel@WrapperModel(model, varargin{:});
            if model.parallel
                % Set parallel properties
                if isempty(model.pool)
                    model.pool = gcp();
                    model.NumWorkers = model.pool.NumWorkers;
                end
                model = validateModel@WrapperModel(model, varargin{:});
                m   = model;
                m.G = [];
                rmodel = getReservoirModel(m);
                if isfield(rmodel, 'FacilityModel')
                    rmodel.FacilityModel = [];
                end
                spmd
                    % Distribute model to all workers
                    lm = m;
                    lm.G = lm.parentModel.G;
                end
                model.localModel = lm;
                model = model.updateSubdomainSetupParallel();
            else
                model = model.updateSubdomainSetupSerial();
            end
        end
        
        %-----------------------------------------------------------------%
        function [model, state] = prepareTimestep(model, varargin)
            % Prepare model and state for timestep
            [model, state] = prepareTimestep@WrapperModel(model, varargin{:});
            % Get partition and update subdomain setups
            p0 = model.partition.value;
            model.partition = model.partition.get(model, varargin{:});
            if isempty(p0) || numel(p0) ~= numel(model.partition.value) ...
                           || any(p0 ~= model.partition.value)
                if ~model.parallel
                    model = model.updateSubdomainSetupSerial();
                else
                    model = model.updateSubdomainSetupParallel();
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function submodel = setSubFacilityModel(model, submodel, mappings) %#ok
            % Set facility model of a sumodel from an already set up
            % facility model for the full model
            rmodel = getReservoirModel(submodel);
            fm = rmodel.FacilityModel;
            fm.ReservoirModel = rmodel;
            fm.ReservoirModel.FacilityModel = [];
            fm.WellModels = fm.WellModels(mappings.wells.keep);
            for i = 1:numel(fm.WellModels)
                fm.WellModels{i}.doUpdatePressureDrop = false;
                fm.WellModels{i}.allowControlSwitching = false;
                fm.WellModels{i}.W.cells = mappings.cells.renum(fm.WellModels{i}.W.cells);
            end
            rmodel.FacilityModel = fm;
            submodel = submodel.setReservoirModel(submodel, rmodel);
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            % TODO: Technically only necessary to compute FacilityFluxProps
            rmodel = getReservoirModel(model);
            if rmodel.FacilityModel.outputFluxes && ~isempty(state.wellSol)
                [~, state] = model.parentModel.getEquations(state0, state, dt, drivingForces, 'resOnly', true);
                state = model.reduceState(state, false);
            end
            [state, report] = model.parentModel.updateAfterConvergence(state0, state, dt, drivingForces);
        end
        
        %-----------------------------------------------------------------%
        function submodel = setSubdomainTolerances(model, submodel)
            % Adjust subdomain tolerances
            tolerances = model.subdomainTol;
            if isempty(tolerances)
                return
            end
            rmodel = getReservoirModel(submodel);
            for i = 1:2:numel(tolerances)
                if isprop(rmodel, tolerances{i})
                    assert(tolerances{i+1} <= 1, 'Subdomain tolerance fraction must be <= 1')
                    rmodel.(tolerances{i}) = rmodel.(tolerances{i}).*tolerances{i+1};
                end
            end
            submodel = submodel.setReservoirModel(submodel, rmodel);
        end
        
        %-----------------------------------------------------------------%
        function fun = getAddStatesFun(model, dummyState) %#ok
            % Get function for adding states
            fun = @(in1, in2) addStatesFun(in1, in2, dummyState);
        end
        
    end
end

%-------------------------------------------------------------------------%
% Helpers

%-------------------------------------------------------------------------%
function state = stripState(state)
    % Strip state to reduce communication overhead
    fields = {'sMax', 'iterations', 'FacilityState', ...
        'eos', 'FacilityState', 'switched', 'switchCount', 'dpRel', 'dpAbs'};
    for f = fields
        if isfield(state, f{1})
            state = rmfield(state, f{1});
        end
    end
end

%-------------------------------------------------------------------------%
function domains = getLocalDomains(numWorkers, p, i)
    % Compute index vector for local domains on each worker
    nDomains = repmat(ceil(max(p)/numWorkers), 1, numWorkers);
    extra    = sum(nDomains) - max(p);
    nDomains(end-extra+1:end) = nDomains(end-extra+1:end) - 1;
    domains = [0, cumsum(nDomains)] + 1;
    domains = domains(i):domains(i+1)-1;
end

%-------------------------------------------------------------------------%
function out = addStatesFun(in1, in2, dummyState)
    % Add up two substates
    map1 = in1.mappings;
    map2 = in2.mappings;
    if isempty(map1) || isempty(map2)
        if isempty(map2)
            out = in1;
        else
            out = in2;
        end
        return
    end
    substate1 = in1.state;
    substate2 = in2.state;
    dummyState = mapState(dummyState, substate1, map1);
    dummyState = mapState(dummyState, substate2, map2);
    map      = mergeMappings(map1, map2);
    substate = getSubState(dummyState, map);
    out      = struct('state', substate, 'mappings', map);
end

%-------------------------------------------------------------------------%
function mappings = mergeMappings(mappings1, mappings2)
    % Merge two subdomain mappings
    if isempty(mappings1)
        mappings1 = mappings2;
    end
    if islogical(mappings1) && islogical(mappings2)
        mappings = mappings1 | mappings2;
        return
    elseif iscell(mappings1) && iscell(mappings2)
        mappings = union(mappings1, mappings2);
        return
    elseif isnumeric(mappings1) && isnumeric(mappings2)
        mappings = mappings1;
        return
    end
    flds = fieldnames(mappings1);
    for i = 1:numel(flds)
        map1 = mappings1.(flds{i});
        map2 = mappings2.(flds{i});
        map  = mergeMappings(map1, map2);
        mappings.(flds{i}) = map;
    end
end

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
