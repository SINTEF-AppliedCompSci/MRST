classdef DomainDecompositionModel < WrapperModel
   
    properties
        partition
        subdomainSetup
        strategy               = 'additive'
        overlap                = 0
        getPartition           = @(varargin) varargin{1}.partition;
        updateFrequency        = 1
        verboseSubmodel        = 0
        subdomainTolerances    = []
        useGlobalPressureRange = false
        parallel               = false
        storeSubModels         = true;
        % Parallel properties
        localModel
        pool    = [];
        domains = [];
        NumWorkers;
    end
    
    methods
        %-----------------------------------------------------------------%
        function model = DomainDecompositionModel(parent, partition, varargin)
            % Base model initialization
            model = model@WrapperModel(parent);
            % Set partition
            dynamicPartition = isa(partition, 'function_handle');
            if dynamicPartition
                model.getPartition = partition;
                partition          = model.getPartition(model);
            end
            model.partition = partition;
            % Merge options
            model = merge_options(model, varargin{:});
            % Control input
            assert(any(strcmpi(model.strategy, {'additive', 'multiplicative'})), ...
                        'Supported strategies are additive and multiplicative');
            if model.parallel
                if strcmpi(model.strategy, 'multiplicative')
                    warning(['Parallel implemetation does not support '       , ...
                             'multiplicative strategy. Changing to additive']);
                end
            end
            if dynamicPartition || max(partition) > 1000
                model.storeSubModels = false;
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
            % Get flux from pressure state if it exists
            if isfield(state, 'statePressure')
                state.flux = state.statePressure.flux;
            end
            if isfield(state, 'FractionalDerivatives')
                state = rmfield(state, 'FractionalDerivatives');
            end
            % Use global pressure range when computing dpRela
            if model.useGlobalPressureRange
                range = max(state.pressure) - min(state.pressure);
                rmodel = model.getReservoirModel();
                if isprop(rmodel, 'incTolPressure')
                    tol   = rmodel.incTolPressure;
                    if isinf(tol)
                        tol = 1e-3;
                    end
                end
                range = max(range, mean(state.pressure)*tol);
                state.pressureRange = range;
            end
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = solveSubDomains(model, state0, dt, drivingForces, state)
            if ~model.parallel
                [state, report] = model.solveSubDomainsSerial(state0, dt, drivingForces, state);
            else
                [state, report] = model.solveSubDomainsParallel(state0, dt, drivingForces, state);
            end
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = solveSubDomainsSerial(model, state0, dt, drivingForces, state)
            % Initialize
            iterations = nan(model.G.cells.num,1);
            [stateInit, stateFinal] = deal(state);
            for i = 1:max(model.partition)
                % Solve subdomain
                [stateInit, stateFinal, iterations] = model.solveSubDomain(model.subdomainSetup{i}, state0, dt, drivingForces, stateInit, stateFinal, iterations);
            end
            state = stateFinal;
            % Make step report
            report = model.makeDomainStepReport('Iterations', iterations);
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = solveSubDomainsParallel(model, state0, dt, drivingForces, state)
            % Initialize
            stateInit = stripState(state);
            state0    = stripState(state0);
            sds = model.subdomainSetup;
            lm  = model.localModel;
            add = @(in1, in2) model.addStates(in1, in2, stateInit);
            ticBytes(model.pool);
            spmd
                localTimer = tic(); % Time the local solve (db only)
                map        = [];
                stateFinal = stateInit;
                iterations = zeros(lm.G.cells.num,1);
                for i = 1:numel(lm.domains)
                    % Solve subdomain
                    [stateInit, stateFinal, iterations, submodel] = lm.solveSubDomain(sds{i}, state0, dt, drivingForces, stateInit, stateFinal, iterations);
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
            report = model.makeDomainStepReport('Iterations', iterations);
        end
        
        
        %-----------------------------------------------------------------%
        function [stateInit, stateFinal, iterations, varargout] = solveSubDomain(model, setup, state0, dt, drivingForces, stateInit, stateFinal, iterations)
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
                    % Multiplicative: update initial state
                    [stateInit, stateFinal] = deal(state);
            end
            % Update iterations
            cells = model.partition == setup.Number;
            iterations(cells) = subreport.Iterations;
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
        function report = makeDomainStepReport(model, varargin)%#ok
            report = struct('Iterations', [], 'WallTime', []);
            report = merge_options(report, varargin{:});
        end
        
        %-----------------------------------------------------------------%
        function model = updateSubdomainSetupSerial(model, partition)
            partition0 = model.partition;
            model.partition = partition;
            np = accumarray(partition,1);
            setup = cell(max(partition),1);
            if all(np == 1) && ~isempty(model.subdomainSetup)
                setup(partition) = model.subdomainSetup(partition0);
            else
                for i = 1:max(partition)
                    setup{i} = model.getSubdomainSetup(i);
                end
            end
            model.subdomainSetup = setup;
        end
        
        %-----------------------------------------------------------------%
        function model = updateSubdomainSetupParallel(model, partition)
            lm = model.localModel;
            spmd
                lm.partition = partition;
                dom = getLocalDomains(lm.NumWorkers, lm.partition, labindex);
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
        function setup = getSubdomainSetup(model, i, construct)
            setup = struct('Model', [], 'NonlinearSolver', [], 'Number', i);
            if nargin < 3
                construct = false;
            end
            if max(model.partition) > 500
                return
            end
            % Make submodel
            verbose  = model.verboseSubmodel;
            cells    = model.partition == i;
            submodel = SubdomainModel(model.parentModel, cells, 'verbose', verbose == 2);
            % Adjust submodel tolerances
            submodel = model.setSubdomainTolerances(submodel);
            % Get linear solver
            rmodel = submodel.getReservoirModel();
            ncomp  = rmodel.getNumberOfComponents();
            ndof   = rmodel.G.cells.num*ncomp;
            if isa(submodel.parentModel, 'ReservoirModel')
                lsol = selectLinearSolverAD(rmodel, 'verbose', verbose);
            elseif isa(submodel.parentModel, 'PressureModel')
                if ndof > 1e4
                    lsol = AMGCLSolverAD('tolerance', 1e-4, 'verbose', verbose);
                else
                    lsol = BackslashSolverAD();
                end
            else
                lsol = BackslashSolverAD();
            end
            isTransport = isa(submodel.parentModel, 'TransportModel');
            % Get nonlinear solver
            nls = NonLinearSolver('minIterations'  , 0                         , ...
                                  'maxIterations'  , 25                        , ...
                                  'ErrorOnFailure' , false                     , ...
                                  'verbose'        , verbose                   , ...
                                  'useRelaxation'  , isTransport               , ...
                                  'LinearSolver'   , lsol                      , ...
                                  'identifier'     , ['SUBDOMAIN ', num2str(i)]);
            % Make subdomain setup
            setup.Model = submodel;
            setup.NonlinearSolver = nls;
        end
        
        %-----------------------------------------------------------------%
        function model = validateModel(model, varargin)
            model = validateModel@WrapperModel(model, varargin{:});
            if model.parallel
                if isempty(model.pool)
                    model.pool = gcp();
                    model.NumWorkers = model.pool.NumWorkers;
                end
                model = validateModel@WrapperModel(model, varargin{:});
                m   = model;
                m.G = [];
                rmodel = m.getReservoirModel();
                if isfield(rmodel, 'FacilityModel')
                    rmodel.FacilityModel = [];
                end
                spmd
                   lm = m;
                   lm.G = lm.parentModel.G;
                end
                model.localModel = lm;
                model = model.updateSubdomainSetupParallel(model.partition);
            else
                model = model.updateSubdomainSetupSerial(model.partition);
            end
        end
        
        %-----------------------------------------------------------------%
        function [model, state] = prepareTimestep(varargin)
            [model, state] = prepareTimestep@WrapperModel(varargin{:});
            % Get partition and update subdomain setups
            part = model.getPartition(varargin{:});
            if any(part ~= model.partition)
                if ~model.parallel
                    model = model.updateSubdomainSetupSerial(part);
                else
                    model = model.updateSubdomainSetupParallel(part);
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function submodel = setSubFacilityModel(model, submodel, mappings) %#ok
            rmodel = submodel.getReservoirModel();
            fm = rmodel.FacilityModel;
            fm.ReservoirModel = rmodel;
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
            [~, state] = model.parentModel.getEquations(state0, state, dt, drivingForces, 'resOnly', true);
            state = model.reduceState(state, false);
            [state, report] = model.parentModel.updateAfterConvergence(state0, state, dt, drivingForces);
        end
        
        %-----------------------------------------------------------------%
        function submodel = setSubdomainTolerances(model, submodel)
            tolerances = model.subdomainTolerances;
            if isempty(tolerances)
                return
            end
            rmodel = submodel.getReservoirModel();
            for i = 1:numel(tolerances.names)
                if isprop(rmodel, tolerances.names{i})
                    rmodel.(tolerances.names{i}) = rmodel.(tolerances.names{i}).*tolerances.factors(i);
                end
            end
            submodel = submodel.setReservoirModel(submodel, rmodel);
        end
        
        %-----------------------------------------------------------------%
        function out = addStates(model, in1, in2, dummyState) %#ok
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
        
    end
end

%-------------------------------------------------------------------------%
% Helpers

%-------------------------------------------------------------------------%
function state = stripState(state)
    fields = {'sMax', 'iterations', 'FacilityState', 'K', ...
        'eos', 'FacilityState', 'switched', 'switchCount', 'dpRel', 'dpAbs'};
    for f = fields
        if isfield(state, f{1})
            state = rmfield(state, f{1});
        end
    end
end

%-------------------------------------------------------------------------%
function domains = getLocalDomains(numWorkers, p, i)
    nDomains = repmat(ceil(max(p)/numWorkers), 1, numWorkers);
    extra    = sum(nDomains) - max(p);
    nDomains(end-extra+1:end) = nDomains(end-extra+1:end) - 1;
    domains = [0, cumsum(nDomains)] + 1;
    domains = domains(i):domains(i+1)-1;
end

%-------------------------------------------------------------------------%
function mappings = mergeMappings(mappings1, mappings2)
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