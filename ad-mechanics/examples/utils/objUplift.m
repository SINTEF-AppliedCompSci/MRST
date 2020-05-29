function obj = objUplift(model, states, schedule, topnode, varargin)
% Compute the average of the vertical displacement at the top of the domain 
% This function is used in runAdjointExample

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('tStep', [], 'ComputePartials', false);
    opt = merge_options(opt, varargin{:});
    
    
    numSteps = numel(schedule.step.val);
    lastStep = numSteps;
    tSteps = opt.tStep;
    if isempty(tSteps) 
        % do all
        tSteps = (1 : numSteps)';
    else
        numSteps = 1;
    end
    obj = repmat({[]}, numSteps, 1);

    for step = 1 : numSteps
        obj{step} = computeUpliftForState(model, states{tSteps(step)}, topnode, ...
                                          'ComputePartials', ...
                                          opt.ComputePartials);
        if tSteps(step) ~= lastStep
            obj{step} = double2ADI(0, obj{step});
        end
    end
    
end
