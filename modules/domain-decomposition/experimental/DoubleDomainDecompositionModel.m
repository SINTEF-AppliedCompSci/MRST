classdef DoubleDomainDecompositionModel < DomainDecompositionModel
    
    properties
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function setup = getSubdomainSetup(model, i, compute)
            setup = struct('Model', [], 'NonlinearSolver', [], 'Number', i);
            p = model.partition.value;
            if nargin < 3
                compute = false;
            end
            if max(p) > 500 && ~compute
                return
            end
            % Make submodel
            verbose  = model.verboseSubmodel;
            cells    = p == i;
            submodel = SubdomainModel(model.parentModel, cells, 'overlap', model.overlap, 'verbose', verbose == 2);
            
            % Adjust submodel tolerances
            submodel = model.setSubdomainTolerances(submodel);
            tp = TopologicalFluxPartition('numBlocks', 20);
            subtol = struct('names', {{'toleranceCNV', 'toleranceMB'}}, 'factors', [1, 1]);
            submodel = DomainDecompositionModel(submodel, tp, 'strategy', 'multiplicative', 'verboseSubmodel', 1, 'subdomainTolerances', subtol);
            % Get linear solver
            rmodel = submodel.getReservoirModel();
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
        function [stateInit, stateFinal, iterations, varargout] = solveSubDomain(model, setup, state0, dt, drivingForces, stateInit, stateFinal, iterations)            % Get submodel
            if isempty(setup.Model)
                setup = model.getSubdomainSetup(setup.Number, true);
            end
            submodel = setup.Model;
            mappings = submodel.parentModel.mappings;
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
            cells = model.partition.value == setup.Number;
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

    end
    
end