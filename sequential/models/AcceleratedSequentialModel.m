classdef AcceleratedSequentialModel < SequentialPressureTransportModel
    % Accelerated sequential model inspired by "Nonlinear acceleration of
    % sequential fully implicit (SFI) method for coupled flow and transport
    % in porous media" by Jiang & Tchelepi, 2019
    properties
        accelerationType % Choice for acceleration
        includePressure % Include pressure updates in acceleration strategy
        m % Algorithm parameter
        w % Default relaxation - algorithm parameter
    end
    
    methods
        function model = AcceleratedSequentialModel(pressureModel, transportModel, varargin)
            model = model@SequentialPressureTransportModel(pressureModel, transportModel);
            model.accelerationType = 'aitken';
            model.m = 3;
            model.w = 0.8;
            model.includePressure = false;
            model.stepFunctionIsLinear = false;
            model.maxOuterIterations = inf;
            model = merge_options(model, varargin{:});
        end
        
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
                if transportReport.Iterations == 0
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

            if ~all(converged)
                state = model.modifyUpdate(state, state0, dt, drivingForces, iteration);
            end
            
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
        
        function state = modifyUpdate(model, state, state0, dt, drivingForces, it)
            s = getTransportVariables(model, state);
            if it == 1
                % Do storage, then return
                s0 = getTransportVariables(model, state0);
                state.transport.s = [s0, s];
                state.transport.r = s - s0;
                state.transport.w = model.w;
                return
            end
            
            isQN = strcmpi(model.accelerationType, 'quasi-newton');
            isAnderson = strcmpi(model.accelerationType, 'anderson');
            isAitken = strcmpi(model.accelerationType, 'aitken');
            isAitkenLocal = strcmpi(model.accelerationType, 'aitken-local');
            
            assert(isQN || isAnderson || isAitken || isAitkenLocal);
            
            s_prev = state.transport.s(:, end);
            r = s - s_prev;
            
            m_v = min(model.m, it);
            lw_bnd = it-m_v+1;
            S = state.transport.s(:, lw_bnd:end);
            R = [state.transport.r(:, lw_bnd:end), r];
            
            if isQN || isAnderson
                dS = diff(S, 1, 2);
                dR = diff(R, 1, 2);
                if isQN
                    gamma = (dS'*dR)\(dS'*r);
                else
                    [Q,R] = qr(dR, 0);
                    gamma = R\(Q'*r);
                end
                gamma(~isfinite(gamma)) = 0;
                s_new = s_prev + model.w*r - (dS + model.w*dR)*gamma;
            elseif isAitken || isAitkenLocal
                wi = state.transport.w;
                r_prev = state.transport.r(:, end);
                if isAitkenLocal
                    w_next = -wi.*r_prev./(r-r_prev);
                else
                    w_next = -wi*(r_prev'*(r - r_prev))/sum((r-r_prev).^2);
                end
                if ~isfinite(w_next)
                    w_next = model.w;
                end
                s_new = s_prev.*(1-w_next) + s.*w_next;
                state.transport.w = w_next;
            end
            
            state = model.setTransportVariables(state, s_new);
            s_new = model.getTransportVariables(state);
            % Storage
            state.transport.s = [state.transport.s, s_new];
            state.transport.r = [state.transport.r, r];
        end
        
        function s = getTransportVariables(model, state)
            if model.includePressure
                s = [state.pressure];
            else
                s = [];
            end
            if isa(model.pressureModel, 'ThreePhaseCompositionalModel')
                s = [s; state.components(:, 1:end-1)];
                if model.water
                    s = [s; state.s(:, 1)];
                end
            else
                s = [s; state.s(:, 1)];
            end
        end
        
        function state = setTransportVariables(model, state, s)
            if model.includePressure
                nc = numel(state.pressure);
                state.pressure = s(1:nc);
                s = s(nc+1:end);
            end
            if isa(model.pressureModel, 'ThreePhaseCompositionalModel')
                n = size(state.components, 2);
                nc = size(state.components, 1);
                
                z = reshape(s(1:(n-1)*nc), [], n-1);
                z = min(max([z, 1 - sum(z, 2)], 0), 1);
                z = z./sum(z, 2);
                state.components = z;
                state = model.pressureModel.computeFlash(state, inf);
                if model.water
                    sw = reshape(s((n-1*nc)+1:n*nc), [], 1);
                    sw = min(max(sw, 0), 1);
                    state.s(:, 1) = sw;
                    state.s = bsxfun(@rdivide, state.s, sum(state.s, 2));
                end
            else
                s = min(max(s, 0), 1);
                state.s(:, 1) = s;
                state.s(:, 2) = 1 - s;
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
