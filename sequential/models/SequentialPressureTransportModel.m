classdef SequentialPressureTransportModel < ReservoirModel
    % Sequential meta-model which solves pressure and transport 
    properties
        pressureModel % Model for computing the pressure
        transportModel % Model for the transport subproblem after pressure is found
        parentModel % Parent class. Used for outer loop convergence if present
        pressureNonLinearSolver % NonLinearSolver instance used for pressure updates
        transportNonLinearSolver % NonLinearSolver instance used for saturation/mole fraction updates

        reupdatePressure % Update pressure based on new mobilities before proceeding to next step
        
        % The following are convergence criteria for the outer loop. To
        % enable outer iterations, set stepFunctionIsLinear = false.
        % Otherwise a single step will be performed. Note that when
        % maxOuterIterations is reached, the solver will act as if the step
        % converged and proceed to the next step.
        maxOuterIterations = inf % Maximum outer loops for a given step before assuming convergence.
        volumeDiscrepancyTolerance = 1e-3; % Tolerance for volume error (|sum of saturations - 1|)
        incTolSaturation = 1e-3; % Tolerance for change in saturations in transport for convergence. 
        outerCheckParentConvergence = true; % Use parentModel, if present, to check convergence
        outerCheckPressureConvergence = false; % Check the pressure equation after transport
        relaxationThreshold = inf; % Outer iteration threshold before we enforce relaxation.
                                   % Setting this to < inf will relpace any relaxation computed by
                                   % the outer nonlinear solver with nls.minRelaxation.
    end
    
    methods
        %-----------------------------------------------------------------%
        function model = SequentialPressureTransportModel(pressureModel, transportModel, varargin)
            model = model@ReservoirModel(transportModel.G);
            model.pressureModel  = pressureModel;
            model.transportModel = transportModel;
            % Default: We do not use outer loop.
            model.stepFunctionIsLinear = true;
            model.reupdatePressure = false;
            [model, extra] = merge_options(model, varargin{:});
            solver_prm = struct('transportLinearSolver', [], ...
                                'pressureLinearSolver', []);
            solver_prm = merge_options(solver_prm, extra{:});
            % Transport model determines the active phases
            if isempty(model.parentModel)
                if isa(model.transportModel, 'WrapperModel')
                    % Wrapper model
                    transportModel = transportModel.parentModel;
                end
                model.water = transportModel.water;
                model.oil   = transportModel.oil;
                model.gas   = transportModel.gas;
            else
                model.water = model.parentModel.water;
                model.oil   = model.parentModel.oil;
                model.gas   = model.parentModel.gas;
            end
            clear transportModel
            if isempty(model.pressureNonLinearSolver)
                model.pressureNonLinearSolver = NonLinearSolver();
            end
            model.pressureNonLinearSolver.identifier = 'PRESSURE';
            
            if isempty(model.transportNonLinearSolver)
                model.transportNonLinearSolver = NonLinearSolver('useRelaxation', true);
            end
            model.transportNonLinearSolver.identifier = 'TRANSPORT';
            
            if ~isempty(solver_prm.pressureLinearSolver)
                model.pressureNonLinearSolver.LinearSolver = ...
                                solver_prm.pressureLinearSolver;
            end
            model.pressureNonLinearSolver.maxTimestepCuts = 0;
            model.pressureNonLinearSolver.errorOnFailure = false;

            if ~isempty(solver_prm.transportLinearSolver)
                model.transportNonLinearSolver.LinearSolver = ...
                                solver_prm.transportLinearSolver;
            end
            
            model.outerCheckParentConvergence = ...
                model.outerCheckParentConvergence && ~isempty(model.parentModel);
            if model.outerCheckParentConvergence
                % Parent model decides outer convergence - cutting only
                % transport does not make sense
                model.transportNonLinearSolver.maxTimestepCuts = 0;
                % Checking saturation increment is not necessary
                model.incTolSaturation = inf;
            end
            
            model.transportNonLinearSolver.errorOnFailure = false;
            model.FacilityModel = []; % Handled by subclasses
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = stepFunction(model, state, state0, dt,...
                                                drivingForces, linsolve, nls,...
                                                iteration, varargin)
            state.iteration = iteration;
            [state, pressure_state, pressureReport, transportReport, pressure_ok, transport_ok, forceArg] =...
                model.solvePressureTransport(state, state0, dt, drivingForces, nls, iteration);
            converged_step = pressure_ok && transport_ok;
            converged = converged_step;
            if converged && ~model.stepFunctionIsLinear
                [converged, values, state] = checkOuterConvergence(model, state, state0, dt, drivingForces, iteration, pressure_state);
                if transportReport.Iterations == 0 && ~model.outerCheckParentConvergence
                    % If the transport did not do anything, we are
                    % effectively converged, even if the values of the
                    % outer residual are not converged. This must be
                    % specifically enabled by allowing zero iterations for
                    % the transport solver and is primarily useful when
                    % there is no reasonable outer convergence criterion.
                    converged = converged | true;
                end
            else
                % Need to have some value in the report
                values = pressureReport.StepReports{end}.NonlinearReport{end}.Residuals(1);
            end
            % If the sequential step failed and we have outer loop enabled,
            % we should leave it up to the outer nonlinear solver what to
            % do.
            failure = ~converged_step && ~model.stepFunctionIsLinear;
            if failure
                FailureMsg = 'Unable to converge transport for given time-step.';
            else
                FailureMsg = '';
            end
            if ~pressure_ok
                converged = converged && false;
            end
            report = model.makeStepReport(...
                                    'Failure',            failure, ...
                                    'Converged',          all(converged), ...
                                    'FailureMsg',         FailureMsg, ...
                                    'ResidualsConverged', converged, ...
                                    'Residuals',          values ...
                                    );
                                
            report.PressureSolver =  pressureReport;
            report.TransportSolver = transportReport;
            if model.reupdatePressure && converged
                state = ...
                    model.pressureNonLinearSolver.solveTimestep(state0, dt, model.pressureModel,...
                            'initialGuess', state, ...
                            forceArg{:});
                [~, state] = model.transportModel.getEquations(state0, state, dt, drivingForces, 'resOnly', true, 'iteration', inf);
                state = model.reduceState(state);
            end
            if isfield(state, 'iteration')
                state = rmfield(state, 'iteration');
            end
        end
        
        %-----------------------------------------------------------------%
        function [state, pressure_state, pressureReport, transportReport, pressure_ok, transport_ok, forceArg] = solvePressureTransport(model, state, state0, dt, drivingForces, nls, iteration)
           % Solve pressure and transport sequentially
            psolver = model.pressureNonLinearSolver;
            tsolver = model.transportNonLinearSolver;
            % Disable connection pressure drop update after first iteration
            prmodel = getReservoirModel(model.pressureModel);
            trmodel = getReservoirModel(model.transportModel);
            if iteration > 1 && isprop(prmodel, 'FacilityModel') ...
                             && ~isempty(prmodel.FacilityModel)
                for i = 1:numel(prmodel.FacilityModel.WellModels)
                    prmodel.FacilityModel.WellModels{i}.doUpdatePressureDrop = false;
                    trmodel.FacilityModel.WellModels{i}.doUpdatePressureDrop = false;
                end
                model.pressureModel  = setReservoirModel(model.pressureModel, prmodel);
                model.transportModel = setReservoirModel(model.transportModel, trmodel);
            end
            % Get the forces used in the step
            forceArg = model.pressureModel.getDrivingForces(drivingForces);
            
            % First, solve the pressure using the pressure nonlinear
            % solver.
            [state, pressureReport] = ...
                psolver.solveTimestep(state0, dt, model.pressureModel,...
                            'initialGuess', state, ...
                            forceArg{:});
            pressure_state = state;
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
                if iteration > model.relaxationThreshold
                    % Enforce relaxation if we are past threshold iteration
                    nls.relaxationParameter = nls.minRelaxation;
                end
                if ~model.stepFunctionIsLinear && ... % Outer loop enabled
                    nls.relaxationParameter ~= 1 % And we should relax
                    % Apply relaxation to global increment if outer solver
                    % requires it
                    state = model.applyRelaxation(state, pressure_state, nls);
                end
                transport_ok = transportReport.Converged;
            else
                transport_ok = false;
                transportReport = [];
            end 
        end
        
        %-----------------------------------------------------------------%
        function state = applyRelaxation(model, state, pressureState, nls)
            % Figure out which variables where updated in the transport step
            [~, names, origin] = model.transportModel.getStateAD(state, false);
            % Exclude everything that is not reservoir variables
            keep = strcmpi(origin, class(getReservoirModel(model.transportModel))) ...
                               | strcmpi(origin, 'TransportModel');    
            names = names(keep);
            % Add pressure if we solver well equations during transport
            if isprop(model.transportModel, 'implicitType') ...
                && strcmpi(model.transportModel.implicitType, 'wells')
                names = [names, 'pressure'];
            end
            % Make sure all satuartions are relaxed
            isS = cellfun(@(name) any(strcmpi(name, {'sW', 'sO', 'sG'})), names);
            if any(isS)
                names = names(~isS);
                names = [names, 's'];
            end
            % Handle blackoil-specific variable switching
            isX = strcmpi(names, 'x');
            if any(isX)
                names = names(~isX);
                names = [names, 'rs', 'rv'];
            end
            % Apply relaxation
            w = nls.relaxationParameter;
            for n = names
                name = n{1};
                vp = model.transportModel.getProp(pressureState, name);
                vt = model.transportModel.getProp(state, name);
                v  = (1-w)*vp + w*vt;
                state = model.transportModel.setProp(state, name, v);
            end
        end
        
        %-----------------------------------------------------------------%
        function [converged, values, state] = checkOuterConvergence(model, state, state0, dt, drivingForces, iteration, pressure_state)
            resnames = {};
            [converged, values] = deal([]);
            % Check volume discrepancy
            tol_vol = model.volumeDiscrepancyTolerance;
            if isfinite(tol_vol)
                if isfield(state, 'sT')
                    sT = state.sT;
                else
                    sT = sum(state.s, 2);
                end
                v = norm(sT - 1, inf);
                values(end+1) = v;
                converged(end+1) = v <= tol_vol;
                resnames{end+1} = 'Volume error';
            end
            % Check increment tolerance
            tol_inc = model.incTolSaturation;
            if isfinite(tol_inc)
                s0 = pressure_state.s;
                s = state.s;
                ds = max(max(abs(s-s0), [], 2));
                values(end+1) = ds;
                converged(end+1) = ds <= tol_inc;
                resnames{end+1} = 'Saturation increment';
            end
            % Alternate mode: If outer loop is enabled, we will revisit
            % the pressue equation to verify that the equation remains
            % converged after the transport step. This check ensures
            % that the assumption of fixed total velocity is reasonable
            % up to some tolerance.
            if ~isempty(model.parentModel) && model.outerCheckParentConvergence
                rmodel = getReservoirModel(model.parentModel);
                if isa(rmodel, 'ThreePhaseCompositionalModel')
                    % Disable pressure increment convergence measure
                    rmodel.incTolPressure = inf;
                    model.parentModel = setReservoirModel(model.parentModel, rmodel);
                end
                state.s = bsxfun(@rdivide, state.s, sum(state.s, 2));
                [problem, state] = model.parentModel.getEquations(state0, state, dt, drivingForces, 'resOnly', true, 'iteration', inf);
                state = model.parentModel.reduceState(state, false);
                [c_parent, v_parent, n_parent] = model.parentModel.checkConvergence(problem);
                converged = [converged, c_parent];
                values = [values, v_parent];
                resnames = [resnames, n_parent];
            end
            
            if model.outerCheckPressureConvergence
                tmp = state;
                tmp.s = bsxfun(@rdivide, state.s, sum(state.s, 2));
                problem = model.pressureModel.getEquations(state0, tmp, dt, drivingForces, 'resOnly', true, 'iteration', inf);
                [c_p, v_p, n_p] = model.pressureModel.checkConvergence(problem);
                keep = ~strcmpi(n_p, 'Delta P');
                converged = [converged, c_p(keep)];
                values = [values, v_p(keep)];
                resnames = [resnames, n_p(keep)];
            end

            converged = converged | iteration > model.maxOuterIterations;
            if model.verbose
                printConvergenceReport(resnames, values, converged, iteration, converged);
            end
        end
        
        %-----------------------------------------------------------------%
        function varargout = getActivePhases(model)
            % Transport model solves for saturations, so that is where the
            % active phases are defined
            varargout = cell(1, nargout);
            [varargout{:}] = model.transportModel.getActivePhases();
        end
        
        %-----------------------------------------------------------------%
        function state = validateState(model, state)
            % Pressure comes first, so validate that.
            state = model.pressureModel.validateState(state);
            state = model.transportModel.validateState(state);
        end

        %-----------------------------------------------------------------%
        function [model, state] = updateForChangedControls(model, state, forces)
            [model.pressureModel, state] = model.pressureModel.updateForChangedControls(state, forces);
            model.transportModel         = model.transportModel.updateForChangedControls(state, forces);
            if ~isempty(model.parentModel)
                model.parentModel = model.parentModel.updateForChangedControls(state, forces);
            end
        end

        %-----------------------------------------------------------------%
        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
            [model.pressureModel, state] = model.pressureModel.prepareTimestep(state, state0, dt, drivingForces);
            model.transportModel         = model.transportModel.prepareTimestep(state, state0, dt, drivingForces);
            if ~isempty(model.parentModel)
                model.parentModel = model.parentModel.prepareTimestep(state, state0, dt, drivingForces);
            end
        end

        %-----------------------------------------------------------------%
        function [model, state] = prepareReportstep(model, state, state0, dt, drivingForces)
            [model.pressureModel, state] = model.pressureModel.prepareReportstep(state, state0, dt, drivingForces);
            model.transportModel         = model.transportModel.prepareReportstep(state, state0, dt, drivingForces);
            if ~isempty(model.parentModel)
                model.parentModel = model.parentModel.prepareReportstep(state, state0, dt, drivingForces);
            end
        end

        %-----------------------------------------------------------------%
        function model = validateModel(model, varargin)
            if isprop(model.pressureModel, 'extraWellSolOutput')
                model.pressureModel.extraWellSolOutput = true;
            end
            model.pressureModel = model.pressureModel.validateModel(varargin{:});
            model.transportModel = model.transportModel.validateModel(varargin{:});
            if ~isempty(model.parentModel)
                model.parentModel = model.parentModel.validateModel(varargin{:});
            end
            return
        end

        %-----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
            [fn, index] = model.pressureModel.getVariableField(name, varargin{:});
        end
        
        %-----------------------------------------------------------------%
        function dt = getMaximumTimestep(model, state, state0, dt, drivingForces)
            dt_p = model.pressureModel.getMaximumTimestep(state, state0, dt, drivingForces);
            dt_t = model.transportModel.getMaximumTimestep(state, state0, dt, drivingForces);
            dt = min(dt_p, dt_t);
        end
        
        %-----------------------------------------------------------------%
        function forces = validateDrivingForces(model, forces, varargin)
           forces = model.pressureModel.validateDrivingForces(forces, varargin{:});
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            if ~isempty(model.parentModel) && ~model.stepFunctionIsLinear ...
                && model.outerCheckParentConvergence
                [state, report] = model.parentModel.updateAfterConvergence(state0, state, dt, drivingForces);
            else
                report = [];
            end
        end
        
    end
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
