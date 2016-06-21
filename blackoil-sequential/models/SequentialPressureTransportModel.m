classdef SequentialPressureTransportModel < ReservoirModel
    properties
        pressureModel
        transportModel
        
        pressureNonLinearSolver
        transportNonLinearSolver
        
        pressureLinearSolver
        transportLinearSolver
        
        outerTolerance
        outerCheckWellConvergence
        maxOuterIterations
    end
    
    methods
        function model = SequentialPressureTransportModel(pressureModel, transportModel, varargin)
            model = model@ReservoirModel([]);
             
            model.pressureModel  = pressureModel;
            model.transportModel = transportModel;
            model.outerTolerance = 1e-3;
            model.outerCheckWellConvergence = false;
            model.maxOuterIterations = 2;
            
            model = merge_options(model, varargin{:});
            
            % Transport model determines the active phases
            model.water = model.transportModel.water;
            model.oil   = model.transportModel.oil;
            model.gas   = model.transportModel.gas;
            
            if isempty(model.pressureNonLinearSolver)
                model.pressureNonLinearSolver = NonLinearSolver();
            end
            model.pressureNonLinearSolver.identifier = 'PRESSURE';
            
            if isempty(model.transportNonLinearSolver)
                model.transportNonLinearSolver = NonLinearSolver();
            end
            model.transportNonLinearSolver.identifier = 'TRANSPORT';
            
            if ~isempty(model.pressureLinearSolver)
                model.pressureNonLinearSolver.LinearSolver = ...
                                model.pressureLinearSolver;
            end
            
            if ~isempty(model.transportLinearSolver)
                model.transportNonLinearSolver.LinearSolver = ...
                                model.transportLinearSolver;
            end
            
            model.transportNonLinearSolver.errorOnFailure = false;
            model.stepFunctionIsLinear = true;
            
        end
        
        function [state, report] = stepFunction(model, state, state0, dt,...
                                                drivingForces, linsolve, nonlinsolve,...
                                                iteration, varargin)
            % Solve pressure
            psolver = model.pressureNonLinearSolver;
            tsolver = model.transportNonLinearSolver;
            
            forceArg = getDrivingForces(model.pressureModel, drivingForces);
            
            [state, pressureReport] = ...
                psolver.solveTimestep(state0, dt, model.pressureModel,...
                            'initialGuess', state, ...
                            forceArg{:});
            pressure_ok = pressureReport.Converged;
            
            if pressure_ok
                % Solve transport
                [state, transportReport] = ...
                    tsolver.solveTimestep(state0, dt, model.transportModel,...
                                'initialGuess', state, ...
                                forceArg{:});
                transport_ok = transportReport.Converged;
            else
                transport_ok = false;
                transportReport = [];
            end
            values = pressureReport.StepReports{end}.NonlinearReport{end}.Residuals;
            
            converged = pressure_ok && transport_ok;
            if converged && ~model.stepFunctionIsLinear
                % Alternate mode: If outer loop is enabled, we will revisit
                % the pressue equation to verify that the equation remains
                % converged after the transport step. This check ensures
                % that the assumption of fixed total velocity is reasonable
                % up to some tolerance.
                problem = model.pressureModel.getEquations(state0, state, dt, drivingForces, 'resOnly', true, 'iteration', inf);
                % Is the pressure still converged when accounting for the
                % updated quantities after transport (mobility, density
                % and so on?)
                [~, values] = model.pressureModel.checkConvergence(problem);
                if model.outerCheckWellConvergence
                    converged = all(values < model.outerTolerance);
                    lv = max(values);
                else
                    converged = values(1) < model.outerTolerance;
                    lv = values(1);
                end
                converged = converged || iteration > model.maxOuterIterations;
                if model.verbose
                    if converged
                        s = 'Converged.';
                    else
                        s = 'Not converged.';
                    end
                    fprintf('OUTER LOOP step #%d with tolerance %1.4e: Largest value %1.4e -> %s \n', ...
                                                            iteration, model.outerTolerance, lv, s);
                end
            end
            if ~pressure_ok
                FailureMsg = 'Pressure failed to converge.';
                failure = true;
            else
                failure = false;
                FailureMsg = '';
            end
            report = model.makeStepReport(...
                                    'Failure',         failure, ...
                                    'Converged',       converged, ...
                                    'FailureMsg',      FailureMsg, ...
                                    'Residuals',       values ...
                                    );
                                
            report.PressureSolver =  pressureReport;
            report.TransportSolver= transportReport;
        end
        
        function varargout = getActivePhases(model)
            varargout = cell(1, nargout);
            [varargout{:}] = model.transportModel.getActivePhases();
        end
        
        function state = validateState(model, state)
            % Pressure comes first, so validate that.
            state = model.pressureModel.validateState(state);
        end
    end
end
