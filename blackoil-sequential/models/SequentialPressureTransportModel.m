classdef SequentialPressureTransportModel < ReservoirModel
    properties
        pressureModel
        transportModel
        
        pressureNonLinearSolver
        transportNonLinearSolver
        
        pressureLinearSolver
        transportLinearSolver
    end
    
    methods
        function model = SequentialPressureTransportModel(pressureModel, transportModel, varargin)
            model = model@ReservoirModel([], [], []);
             
            model.pressureModel  = pressureModel;
            model.transportModel = transportModel;
            
            model = merge_options(model, varargin{:});
            
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
            
            model.stepFunctionIsLinear = true;
        end
        
        function [state, report] = stepFunction(model, state, state0, dt,...
                                                drivingForces, linsolve, nonlinsolve,...
                                                iteration, varargin)
            % Solve pressure
            psolver = model.pressureNonLinearSolver;
            tsolver = model.transportNonLinearSolver;
            
            [state, pressureReport] = ...
                psolver.solveTimestep(state0, dt, model.pressureModel,...
                            'Wells', drivingForces.Wells, ...
                            'bc', drivingForces.bc, ...
                            'src', drivingForces.src);
            pressure_ok = pressureReport.Converged;
            
            if pressure_ok
                % Solve transport
                [state, transportReport] = ...
                    tsolver.solveTimestep(state, dt, model.transportModel,...
                                'Wells', drivingForces.Wells, ...
                                'bc', drivingForces.bc, ...
                                'src', drivingForces.src);
                transport_ok = transportReport.Converged;
            else
                transport_ok = false;
                transportReport = [];
            end
            converged = pressure_ok && transport_ok;
            
            if ~pressure_ok
                FailureMsg = 'Pressure failed to converge!';
            elseif ~transport_ok
                FailureMsg = 'Transport failed to converge!';
            else
                FailureMsg = '';
            end

            values = [];
            report = model.makeStepReport(...
                                    'Failure',         ~converged, ...
                                    'Converged',       converged, ...
                                    'FailureMsg',      FailureMsg, ...
                                    'Residuals',       values ...
                                    );
                                
            report.PressureSolver =  pressureReport;
            report.TransportSolver= transportReport;
        end
        
        function varargout = getActivePhases(model)
            [varargout{:}] = model.transportModel.getActivePhases(model);
        end
        
    end
end
