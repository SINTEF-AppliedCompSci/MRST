function obj = objUplift(model, states, schedule, topnode, varargin)
% Compute the average of the vertical displacement at the top of the domain 
% This function is used in runAdjointExample

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

    opt = struct('computeAllSteps'      , true , ...
                 'tStep'                , []   , ...
                 'state'                , [], ...
                 'ComputePartials'      , false);
    
    opt = merge_options(opt, varargin{:});
    
    if opt.computeAllSteps
        % Do all the steps
        numSteps      = numel(schedule.step.val);
        tSteps        = (1 : numSteps)';
        scheduleSteps = tSteps;
    else
        % Compute for the given step. When 'computePartial' is used, we need state as input as it contained the AD instantiation.
        numSteps = 1;
        states{1} = opt.state;
        scheduleSteps = opt.tStep;
        tSteps = 1;
    end


    obj = repmat({[]}, numSteps, 1);

    for step = 1 : numSteps
        dt = schedule.step.val(scheduleSteps(step));
        obj{step} = computeUpliftForState(model, states{tSteps(step)}, topnode, 'ComputePartials', ...
                                          opt.ComputePartials);
        obj{step} = obj{step}*dt;
    end
    
end
