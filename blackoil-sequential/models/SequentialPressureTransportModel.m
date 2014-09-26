classdef SequentialPressureTransportModel < ReservoirModel
    properties
        pressureModel
        transportModel
    end
    
    methods
        function model = SequentialPressureTransportModel(pressureModel, transportModel)
            model = model@ReservoirModel([], [], []);
             
            model.pressureModel  = pressureModel;
            model.transportModel = transportModel;
        end
        
        function [state, report] = stepFunction(model, state, state0, dt,...
                                                drivingForces, linsolve, nonlinsolve,...
                                                onlyCheckConvergence, varargin)
            % Solve pressure
            [state, pressureReport] = ...
                nonlinsolve.solveTimestep(state0, dt, model.pressureModel,...
                            'Wells', drivingForces.Wells, ...
                            'bc', drivingForces.bc, ...
                            'src', drivingForces.src);
            pressure_ok = pressureReport.Converged;
            
            if pressure_ok
                % Solve transport
                [state, transportReport] = ...
                    nonlinsolve.solveTimestep(state, dt, model.transportModel,...
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
            
%             if converged
%                 values = transportReport.Residuals;
%             else
                values = [];
%             end
            report = model.makeStepReport(...
                                    'Failure',         ~converged, ...
                                    'Converged',       converged, ...
                                    'FailureMsg',      FailureMsg, ...
                                    'Residuals',       values ...
                                    );
                                
            report.PressureSolver =  pressureReport;
            report.TransportSolver= transportReport;
        end
    end
end