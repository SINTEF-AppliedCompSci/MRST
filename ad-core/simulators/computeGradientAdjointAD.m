function gradients = computeGradientAdjointAD(state0, states, model, schedule, getObjective, varargin)
%Compute gradients using an adjoint/backward simulation that is linear in each step
%
% SYNOPSIS:
%   grad = computeGradientAdjointAD(state0, states, model, schedule, getObjective
%
% DESCRIPTION:
%   For a given schedule, compute gradients of objective with respect to well 
%   controls by solving adjoint equations.
%
% REQUIRED PARAMETERS:
%
%   state0       - Physical model state at t = 0
%
%   states       - All previous states. Must support the syntax 
%                  `state = states{i}`. If the problem is too large to fit in
%                  memory, it can use ResultHandler class to retrieve files
%                  from the disk.
%                  
%   model        - Subclass of `PhysicalModel` class such as
%                  `ThreePhaseBlackOilModel` that models the physical 
%                  effects we want to study.
%
%   schedule     - Schedule suitable for `simulateScheduleAD`.
%
%   getObjective - Function handle for getting objective function value 
%                  for a given timestep with derivatives. Format: @(tstep)
%
% OPTIONAL PARAMETERS:
%
%  'Verbose'        - Indicate if extra output is to be printed such as
%                     detailed convergence reports and so on. 
% 
%  'scaling'        - Struct with fields `rate` and `pressure` used to
%                     scale the relevant control equations, if the model 
%                     supports it.
%
%  'LinearSolver'   - Subclass of `LinearSolverAD` suitable for solving the
%                     adjoint systems.
%
% RETURNS:
%   gradients - Cell array of gradients for each control step.
%
% SEE ALSO:
%   `computeGradientPerturbationAD`, `simulateScheduleAD`

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

    opt = struct('ControlVariables', {{'well'}}, ...
                 'LinearSolver',     [], ...
                 'OutputPerTimestep', false, ...
                 'Verbose', mrstVerbose());
    opt = merge_options(opt, varargin{:});
    
    if isempty(opt.LinearSolver)
        linsolve = BackslashSolverAD();
    else
        linsolve = opt.LinearSolver;
    end
    
    if iscell(opt.ControlVariables)
        ncv = numel(opt.ControlVariables);
    else
        ncv = 1;
    end
    dispif(opt.Verbose, 'Preparing model for simulation...\n')
    ctrl = schedule.control(schedule.step.control(1));
    [~, fstruct] = model.getDrivingForces(ctrl);
    model = model.validateModel(fstruct);
    dispif(opt.Verbose, 'Model ready for simulation...\n')
    
    % Check if intiial state is reasonable
    dispif(opt.Verbose, 'Validating initial state...\n')
    state0 = model.validateState(state0);
    dispif(opt.Verbose, 'Initial state ready for simulation.\n')
    
    % Function handle to pick initial state
    getState = @(i) getStateFromInput(schedule, states, state0, i);
    
    nstep = numel(schedule.step.val);
    grad = [];
    gradstep = cell(nstep, ncv);
    nt = nstep;
    for step = nt:-1:1
        fprintf('Solving reverse mode step %d of %d\n', nt - step + 1, nt);
        [dg, grad, report] = model.solveAdjoint(linsolve, getState, ...
                                         getObjective, schedule, grad, step);
        if isa(model.FacilityModel, 'GenericFacilityModel')
            ws = states{step}.wellSol;
            cntrScale = getControlEqScaling({ws.type}, model.FacilityModel);
            eqNo      = strcmp('well', report.Types);
            st        = vertcat(ws.status);
            dg{eqNo}  = cntrScale(st).*dg{eqNo};
        end
        gradstep(step, :) = getRequestedGradients(dg, report, opt.ControlVariables);
    end
    
    % Sum up to the control steps
    if ~opt.OutputPerTimestep
        nOut = numel(schedule.control);
        contr = schedule.step.control;
    else % output gradient per time-step
        nOut = nt;
        contr = (1:nt)';
    end
    gradients = cell(ncv, nOut);
    
    for k = 1:nOut
        ck = contr == k;
        for j = 1:ncv
            tmp = gradstep(ck, j);
            gradients{j, k} = 0;
            for i = 1:numel(tmp)
                gradients{j, k} = gradients{j, k} + tmp{i};
            end
        end
    end
end

function state = getStateFromInput(schedule, states, state0, i)
    if i == 0
        state = state0;
    elseif i > numel(schedule.step.val)
        state = [];
    else
        state = states{i};
    end
end

function g = getRequestedGradients(dg, report, wantGradFor)
    if ischar(wantGradFor)
        g = {vertcat(dg{strcmpi(report.Types, wantGradFor)})};
    else
        ng = numel(wantGradFor);
        g = cell(1, ng);
        for i = 1:ng
            n = wantGradFor{i};
            g{i} = vertcat(dg{strcmpi(report.Types, n)});
        end
    end
end

function sc = getControlEqScaling(tp, facility)
% In order to produce correct gradients, this scaling must be precisely the
% same as in the control equation. This is rather shaky!
sc    = ones(numel(tp), 1);
is_bhp = strcmp('bhp', tp);
if any(is_bhp)
    if isfinite(facility.toleranceWellRate) && isfinite(facility.toleranceWellBHP)
        eqScale = facility.toleranceWellRate/facility.toleranceWellBHP;
    else
        eqScale = 1/(day()*barsa());
    end
    sc(is_bhp) = eqScale;
end
end

    
