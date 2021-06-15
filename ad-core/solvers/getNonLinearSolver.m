function nonlinear = getNonLinearSolver(model, varargin)
% Set up reasonable defaults for the nonlinear solver for a field
% simulation with significant size and complexity
%
% SYNOPSIS:
%   solver = getNonLinearSolver(model)
%
% REQUIRED PARAMETERS:
%   model  - Simulation model (subclass of PhysicalModel).
%
% OPTIONAL PARAMETERS:
%   'DynamicTimesteps'  - Set up a simple iteration count timestep
%                         selector.
%
%   'useCPR'            - Set up CPR-type preconditioner as the linear
%                         solver. Will try to use the best known linear
%                         solver (either AGMG or Matlab builtin at the
%                         moment).
% RETURNS:
%   solver    - NonLinearSolver instance.
%
% SEE ALSO:
%   simulateScheduleAD, NonLinearSolver

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

    opt = struct('TimestepStrategy',      'iteration', ...
                 'useCPR',                true, ...
                 'LinearSolverArguments', {{}});

    [opt, varg] = merge_options(opt, varargin{:});

    nonlinear = NonLinearSolver(varg{:});
    nonlinear.LinearSolver = selectLinearSolverAD(model, 'useCPR', opt.useCPR, ...
                                                    opt.LinearSolverArguments{:});
    switch lower(opt.TimestepStrategy)
        case 'none'
            % Do nothing
            sel = [];
        case 'iteration'
            % Control on iterations
            sel = IterationCountTimeStepSelector('targetIterationCount', 8);
        case 'ds'
            % Control on saturation change
            sel = ...
                StateChangeTimeStepSelector('targetProps', {'s'},...
                                            'targetChangeAbs', 0.2, ...
                                            'targetIterationCount', inf);
        case 'dsdc'
            % Control on saturation + components
            names = {'s'};
            targets = 0.2;
            if isa(model, 'ThreePhaseCompositionalModel')
                names = {'s', 'components'};
                targets = [0.2, 0.2];
            end
            sel = ...
                StateChangeTimeStepSelector('targetProps', names,...
                                            'targetChangeAbs', targets, ...
                                            'targetIterationCount', inf);

        otherwise
            error('Unknown timestepping strategy %s', opt.TimestepStrategy);
    end
    if ~isempty(sel)
        sel.firstRampupStepRelative = 0.1;
        sel.firstRampupStep = 1*day;
        nonlinear.timeStepSelector = sel;
    end
end
