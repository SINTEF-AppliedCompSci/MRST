function grad = computeGradientPerturbationAD(state0, model, schedule, getObjective, varargin)
%Compute gradients using finite difference approximation by perturbing controls
%
% SYNOPSIS:
%   grad = computeGradientPerturbationAD(state0, model, schedule, getObjective)
%
% DESCRIPTION:
%   For a given schedule, compute gradients with regards to well controls
%   by perturbing all controls ever so slightly and re-running the
%   simulation.
%  
%   As the cost of this routine grows is approximately ::
%
%        (# wells)x(# ctrl step) x cost of schedule
%
%   it can be potentially extremely expensive. It is better to use the
%   `computeGradientAdjointAD` routine for most practical purposes. This
%   routine is primarily designed for validation of said routine.
%
% REQUIRED PARAMETERS:
%
%   state0       - Physical model state at `t = 0`
%
%   model        - Subclass of `PhysicalModel` class such as
%                  `ThreePhaseBlackOilModel` that models the physical 
%                  effects we want to study.
%
%   schedule     - Schedule suitable for `simulateScheduleAD`.
%
%   getObjective - Function handle for getting objective function value 
%                  from a set of wellSols. 
%                  Function handle format: `@(wellSols, states, schedule)`
%
% OPTIONAL PARAMETERS:
%   
%  'Verbose'        - Indicate if extra output is to be printed such as
%                     detailed convergence reports and so on. 
% 
%  'scaling'        - Struct with fields `rate` and `pressure` used to
%                     scale the numerical perturbation.
%
%  'perturbation'   - Magnitude of perturbation. Default `1e-7`.
%
%
% RETURNS:
%   grad - Gradients for each step.
%
% SEE ALSO:
%   `computeGradientAdjointAD`, `simulateScheduleAD`

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

    opt = struct('Verbose', mrstVerbose(),...
                 'perturbation', 1e-7, ...
                 'scaling', []);
    opt = merge_options(opt, varargin{:});

    if ~isempty(opt.scaling)
        scalFacs = opt.scaling;
    else
        scalFacs.rate = 1; scalFacs.pressure = 1;
    end
    
    solve = @(schedule) simulateScheduleAD(state0, model, schedule, 'Verbose', opt.Verbose);
    schedule0 = schedule;

    % Solve entire schedule to get baseline of the objective function
    [wellSols,states] = solve(schedule0);

    % Set up objective function storage
    computeObj = @(ws,states,schedule) sum(cell2mat(getObjective(ws,states,schedule)));
    val0 = computeObj(wellSols,states,schedule0);

    grad = cell(1, numel(schedule0.control));
    % Run a schedule per well, for each control step, perturbing the
    % control variable ever so slightly to get a local finite difference
    % approximation of the gradient.
    for cn = 1:numel(schedule.control)
        
        ctrl = schedule0.control(cn);
        nWell = numel(ctrl.W);
        grad{cn} = zeros(nWell, 1);

        dispif(opt.Verbose, 'Solving for control %d of %d', cn, numel(schedule.control));
        for k = 1:nWell
            dispif(opt.Verbose, 'Processing well %d of %d\n', k, numel(ctrl.W));

            w = ctrl.W(k);

            if strcmpi(w.type, 'bhp')
                e = scalFacs.pressure*opt.perturbation;
            else
                e = scalFacs.rate*opt.perturbation;
            end

            schedule = schedule0;
            w.val = w.val + e;
            schedule.control(cn).W(k) = w;
            
            [wellSols,states] = solve(schedule);
            valk = sum( cell2mat(getObjective(wellSols,states,schedule)));
            grad{cn}(k) = (valk-val0)/e;
        end
    end
end


