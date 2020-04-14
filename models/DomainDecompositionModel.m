classdef DomainDecompositionModel < WrapperModel
   
    properties
        partition
        subdomainSetup
        type            = 'additive'
        overlap         = 0
        getPartition    = @(model, varargin) model.partition
        verboseSubmodel = false
        subdomainTolerances = []
        useGlobalPressureRange = false
    end
    
    methods
        %-----------------------------------------------------------------%
        function model = DomainDecompositionModel(parent, partition, varargin)
            % Base model initialization
            model = model@WrapperModel(parent);
            % Set partition
            model.partition = partition;
            % Merge options
            model = merge_options(model, varargin{:});
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = stepFunction(model, state, state0, dt, drivingForces, linsolver, nonlinsolver, iteration, varargin)
            % Prepare state
            state = model.prepareState(state, iteration);
            % Get residual
            outOfIterations = iteration > nonlinsolver.maxIterations;
            problem = model.getEquations(state0, state, dt, drivingForces, ...
                                       'ResOnly', true, ...
                                       'iteration', iteration, ...
                                       varargin{:});
            problem.iterationNo   = iteration;
            problem.drivingForces = drivingForces;
            % Check convergence
            [convergence, values, resnames] = model.parentModel.checkConvergence(problem);
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
                state.iterations = state.iterations + subdomainReport.Iterations;
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
        function state = prepareState(model, state, iteration)
            % Initialize total subdomain iterations
            if iteration == 1
                state.iterations = zeros(model.G.cells.num,1);
            end
            % Get flux from pressure state if it exists
            if isfield(state, 'statePressure')
                state.flux = state.statePressure.flux;
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
        function [state, report, substates] = solveSubDomains(model, state0, dt, drivingForces, state)
            part = model.partition;
            % Initialize
            iterations = nan(model.G.cells.num,1);
            [initState, finalState] = deal(state);
            substates = cell(max(part),1);
            for i = 1:max(part)
                % Get cell subset
                cells = part == i;
                % Solve subdomain
                [state, substates{i}, subreport, mappings] = model.solveSubDomain(model.subdomainSetup{i}, state0, dt, drivingForces, initState);
                % Update initial and final state
                switch model.type
                    case 'additive'
                        % Additive method: map substate to final state
                        if subreport.Converged
                            finalState = mapState(finalState, substates{i}, mappings);
                        end
                    case 'multiplicative'
                        % Multiplicative: update initial state
                        [initState, finalState] = deal(state);
                end
                % Update iterations
                iterations(cells) = subreport.Iterations;
            end
            state = finalState;
            % Make step report
            report = model.makeDomainStepReport('Iterations', iterations);
        end
        
        %-----------------------------------------------------------------%
        function [state, substate, subreport, mappings, substate0, subforces, submodel] = solveSubDomain(model, setup, state0, dt, drivingForces, state)
            % Get submodel
            submodel = setup.Model;
            mappings = submodel.mappings;
            % Get subdomain forces
            [subforces, mappings] = getSubForces(drivingForces, mappings);
            % Set facility model
            submodel = model.setSubFacilityModel(submodel, mappings);
            % Get subdomain state
            [substateInit, mappings] = getSubState(state, mappings);
            substate0                = getSubState(state0, mappings);
            % Solve timestep
            nls      = setup.NonlinearSolver;
            forceArg = getDrivingForces(submodel, subforces);
            [substate, subreport] = nls.solveTimestep(substate0, dt, submodel, 'initialGuess', substateInit, forceArg{:});
            if ~subreport.Converged
               substate = substateInit; 
            end
            % Map substate to state
            state = mapState(state, substate, mappings);
        end
        
        %-----------------------------------------------------------------%
        function report = makeDomainStepReport(model, varargin)%#ok
            report = struct('Iterations', [], 'WallTime', []);
            report = merge_options(report, varargin{:});
        end
        
        %-----------------------------------------------------------------%
        function model = updateSubdomainSetup(model, partition)
            setup = cell(max(partition),1);
            model.partition = partition;
            for i = 1:max(partition)
                cells = partition == i;
                setup{i} = model.getSubdomainSetup(cells);
            end
            model.subdomainSetup = setup;
        end
        
        %-----------------------------------------------------------------%
        function setup = getSubdomainSetup(model, cells)
            % Make submodel
            submodel = SubdomainModel(model.parentModel, cells, 'verbose', model.verboseSubmodel);
            % Adjust submodel tolerances
            submodel = model.setSubdomainTolerances(submodel);
            % Get linear solver
            rmodel = submodel.getReservoirModel();
            ncomp  = rmodel.getNumberOfComponents();
            ndof   = rmodel.G.cells.num*ncomp;
            if isa(submodel.parentModel, 'ReservoirModel')
                lsol = selectLinearSolverAD(rmodel, 'verbose', model.verboseSubmodel);
            elseif isa(submodel.parentModel, 'PressureModel')
                if ndof > 1e4
                    lsol = AMGCLSolverAD('tolerance', 1e-4, 'verbose', model.verboseSubmodel);
                else
                    lsol = BackslashSolverAD();
                end
            else
                lsol = BackslashSolverAD();
            end
            isTransport = isa(submodel.parentModel, 'TransportModel');
            % Get nonlinear solver
            nls = NonLinearSolver('minIterations'  , 0          , ...
                                  'maxIterations'  , 25         , ...
                                  'ErrorOnFailure' , false      , ...
                                  'verbose'        , false      , ...
                                  'useRelaxation'  , isTransport, ...
                                  'LinearSolver'   , lsol       );
            % Make subdomain setup
            setup = struct('Model', submodel, 'NonlinearSolver', nls);
        end
        
        %-----------------------------------------------------------------%
        function model = validateModel(model, varargin)
            model = validateModel@WrapperModel(model, varargin{:});
            model = model.updateSubdomainSetup(model.partition);
        end
        
        %-----------------------------------------------------------------%
        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
            [model, state] = prepareTimestep@WrapperModel(model, state, state0, dt, drivingForces);
            % Get partition and update subdomain setups
            part = model.getPartition(model, state, state0, dt, drivingForces);
            if any(part ~= model.partition)
                model = model.updateSubdomainSetup(part);
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
        
    end
end