classdef SequentialPressureTransportModel < ReservoirModel
    % Sequential meta-model which solves pressure and transport using a fixed
    % flux splitting
    properties
        % Model for computing the pressure
        pressureModel
        % Model for the transport subproblem after pressure is found
        transportModel
        parentModel
        
        % NonLinearSolver instance used for pressure updates
        pressureNonLinearSolver
        % NonLinearSolver instance used for saturation/mole fraction updates
        transportNonLinearSolver
        % Utility prop for setting pressure linear solver
        pressureLinearSolver
        % Utility prop for setting transport linear solver
        transportLinearSolver
        
        % Outer tolerance, which, if stepFunctionIsLinear is set to false,
        % is used to check if the pressure must be recomputed after
        % transport has been solved, in order to converge to the fully
        % implicit solution.
        outerTolerance
        % Indicates if we check well values when outer loop is enabled
        outerCheckWellConvergence
        % Maximum outer loops for a given step. When maxOuterIterations is
        % reached, the solver will act as if the step converged and
        % continue.
        maxOuterIterations
        % Update pressure based on new mobilities before proceeding to next
        % step
        reupdatePressure
    end
    
    methods
        function model = SequentialPressureTransportModel(pressureModel, transportModel, varargin)
            model = model@ReservoirModel(transportModel.G);
            % Set up defaults
            model.pressureModel  = pressureModel;
            model.transportModel = transportModel;
            model.outerTolerance = 1e-3;
            model.outerCheckWellConvergence = false;
            model.maxOuterIterations = 2;
            % Default: We do not use outer loop.
            model.stepFunctionIsLinear = true;
            model.reupdatePressure = false;
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
            model.pressureNonLinearSolver.maxTimestepCuts = 0;
            model.pressureNonLinearSolver.errorOnFailure = false;

            if ~isempty(model.transportLinearSolver)
                model.transportNonLinearSolver.LinearSolver = ...
                                model.transportLinearSolver;
            end
            
            model.transportNonLinearSolver.errorOnFailure = false;
        end
        
        function [state, report] = stepFunction(model, state, state0, dt,...
                                                drivingForces, linsolve, nonlinsolve,...
                                                iteration, varargin)
            
            [state, pressureReport, transportReport, pressure_ok, transport_ok, forceArg] =...
                model.solvePressureTransport(state, state0, dt, drivingForces, iteration);

            converged = pressure_ok && transport_ok;
            if converged && ~model.stepFunctionIsLinear
                [converged, values, state] = checkOuterConvergence(model, state, state0, dt, drivingForces, iteration);
            else
                % Need to have some value in the report
                values = pressureReport.StepReports{end}.NonlinearReport{end}.Residuals(1);
            end
            failure = false;
            FailureMsg = '';
            if ~pressure_ok
                converged = converged && false;
            end
            report = model.makeStepReport(...
                                    'Failure',         failure, ...
                                    'Converged',       all(converged), ...
                                    'FailureMsg',      FailureMsg, ...
                                    'ResidualsConverged', converged, ...
                                    'Residuals',       values ...
                                    );
                                
            report.PressureSolver =  pressureReport;
            report.TransportSolver = transportReport;
            
            if model.reupdatePressure && converged
                state = ...
                    psolver.solveTimestep(state0, dt, model.pressureModel,...
                            'initialGuess', state, ...
                            forceArg{:});
                [~, state] = model.transportModel.getEquations(state0, state, dt, drivingForces, 'resOnly', true, 'iteration', inf);
            end
        end
        
        function [state, pressureReport, transportReport, pressure_ok, transport_ok, forceArg] = solvePressureTransport(model, state, state0, dt, drivingForces, iteration)
           % Solve pressure and transport sequentially
            psolver = model.pressureNonLinearSolver;
            tsolver = model.transportNonLinearSolver;
            if iteration > 1
                for i = 1:numel(model.pressureModel.FacilityModel.WellModels)
                    model.pressureModel.FacilityModel.WellModels{i}.doUpdatePressureDrop = false;
                end
            end
            % Get the forces used in the step
            forceArg = model.pressureModel.getDrivingForces(drivingForces);
            
            % First, solve the pressure using the pressure nonlinear
            % solver.
            [state, pressureReport] = ...
                psolver.solveTimestep(state0, dt, model.pressureModel,...
                            'initialGuess', state, ...
                            forceArg{:});
            pressure_ok = pressureReport.Converged || psolver.continueOnFailure;
            
            if pressure_ok
                if ~isempty(drivingForces.bc)
                    isDir = strcmpi(drivingForces.bc.type, 'pressure');
                    if any(isDir)
                        % Convert Dirichlet boundary conditions to flux
                        % boundary conditions for the transport
                        transportForces = drivingForces;
                        G = model.pressureModel.G;
                        dirFace = transportForces.bc.face(isDir);
                        q = sum(state.flux(dirFace, :), 2);
                        sgn = 1 - 2*(G.faces.neighbors(dirFace, 2) == 0);
                        transportForces.bc.value(isDir) = sgn.*q;
                        [transportForces.bc.type{isDir}] = deal('resflux');
                        forceArg = model.transportModel.getDrivingForces(transportForces);
                    end
                end
                state.timestep = dt;
                state.pressure_full = state.pressure;
                % If pressure converged, we proceed to solve the transport
                [state, transportReport] = ...
                    tsolver.solveTimestep(state0, dt, model.transportModel,...
                                'initialGuess', state, ...
                                forceArg{:});
                transport_ok = transportReport.Converged;
            else
                transport_ok = false;
                transportReport = [];
            end 
        end
        
        function [converged, values, state] = checkOuterConvergence(model, state, state0, dt, drivingForces, iteration)
            % Alternate mode: If outer loop is enabled, we will revisit
            % the pressue equation to verify that the equation remains
            % converged after the transport step. This check ensures
            % that the assumption of fixed total velocity is reasonable
            % up to some tolerance.
            if ~isempty(model.parentModel)
                state.s = bsxfun(@rdivide, state.s, sum(state.s, 2));
                [problem, state] = model.parentModel.getEquations(state0, state, dt, drivingForces, 'resOnly', true, 'iteration', inf);
                [converged, values, resnames] = model.parentModel.checkConvergence(problem);
                if model.verbose
                    printConvergenceReport(resnames, values, converged, iteration);
                end
            else
                if isa(model.pressureModel, 'ThreePhaseCompositionalModel')
                    values = max(abs(sum(state.s, 2) - 1));
                    if ~model.transportModel.useIncTolComposition
                        % Make a normalization of saturations and check if
                        % the equations are still converged.
                        state_normalized = state;
                        state_normalized.s = bsxfun(@rdivide, state_normalized.s, sum(state_normalized.s, 2));
                        
                        [problem, state_normalized] = model.transportModel.getEquations(state0, state_normalized, dt, drivingForces, 'resOnly', true, 'iteration', inf);
                        conv_t = model.transportModel.checkConvergence(problem);
                        if all(conv_t)
                            state = state_normalized;
                            values = 0;
                        end
                    end
                else
                    problem = model.pressureModel.getEquations(state0, state, dt, drivingForces, 'resOnly', true, 'iteration', inf);
                    % Is the pressure still converged when accounting for the
                    % updated quantities after transport (mobility, density
                    % and so on?)
                    [~, values] = model.pressureModel.checkConvergence(problem);
                end
                if model.outerCheckWellConvergence
                    lv = max(values);
                else
                    values = values(1);
                    lv = values(1);
                end
                converged = all(values < model.outerTolerance);
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
        end
        
        function varargout = getActivePhases(model)
            % Transport model solves for saturations, so that is where the
            % active phases are defined
            varargout = cell(1, nargout);
            [varargout{:}] = model.transportModel.getActivePhases();
        end
        
        function state = validateState(model, state)
            % Pressure comes first, so validate that.
            state = model.pressureModel.validateState(state);
        end

        function [model, state] = updateForChangedControls(model, state, forces)
            [model.pressureModel, state] = model.pressureModel.updateForChangedControls(state, forces);
        end

        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
            [model.pressureModel, state] = model.pressureModel.prepareTimestep(state, state0, dt, drivingForces);
        end

        function model = validateModel(model, varargin)
            model.pressureModel.extraWellSolOutput = true;
            model.pressureModel = model.pressureModel.validateModel(varargin{:});
            model.transportModel = model.transportModel.validateModel(varargin{:});
            if ~isempty(model.parentModel)
                model.parentModel = model.parentModel.validateModel(varargin{:});
            end
            return
        end

        function [fn, index] = getVariableField(model, name)
            [fn, index] = model.pressureModel.getVariableField(name);
        end
    end
end

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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
