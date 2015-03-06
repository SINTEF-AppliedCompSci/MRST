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
            
            if isempty(model.transportNonLinearSolver)
                model.transportNonLinearSolver = NonLinearSolver();
            end
            
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
            
            fprintf('@@@ PRESSURE STEP\n');
            
            [state, pressureReport] = ...
                psolver.solveTimestep(state0, dt, model.pressureModel,...
                            'initialGuess', state, ...
                            'Wells', drivingForces.Wells, ...
                            'bc', drivingForces.bc, ...
                            'src', drivingForces.src);
            pressure_ok = pressureReport.Converged;
            
            if pressure_ok
                
                fprintf('@@@ TRANSPORT STEP\n');
                
                % Solve transport
                [state, transportReport] = ...
                    tsolver.solveTimestep(state0, dt, model.transportModel,...
                                'initialGuess', state, ...
                                'Wells', drivingForces.Wells, ...
                                'bc', drivingForces.bc, ...
                                'src', drivingForces.src);
                transport_ok = transportReport.Converged;
            else
                
                fprintf('@@@ pressure step not ok\n');
                
                transport_ok = false;
                transportReport = [];
            end
            values = pressureReport.StepReports{end}.NonlinearReport{end}.Residuals;
            
            converged = pressure_ok && transport_ok;
            if converged && ~model.stepFunctionIsLinear
                problem = model.pressureModel.getEquations(state0, state, dt, drivingForces, 'resOnly', true, 'iteration', inf);
                % Is the pressure still converged when accounting for the
                % mobilities?
                [~, values] = model.pressureModel.checkConvergence(problem);
                if model.outerCheckWellConvergence
                    converged = all(values < model.outerTolerance);
                else
                    converged = values(1) < model.outerTolerance;
                end
                converged = converged || iteration > model.maxOuterIterations;
            end
            if ~pressure_ok
                FailureMsg = 'Pressure failed to converge!';
            else
                FailureMsg = '';
            end

            
            report = model.makeStepReport(...
                                    'Failure',        ~pressure_ok, ...
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
        
    end
end
