function solver = getNonLinearSolver(model, varargin)
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

opt = struct('DynamicTimesteps', true, ...
             'useCPR',           true);

[opt, varg] = merge_options(opt, varargin{:});

[linsolve, timestepper] = deal([]);
if opt.useCPR && isa(model, 'ReservoirModel')
    if ~isempty(mrstPath('agmg'))
        mrstModule add agmg
        pSolver = AGMGSolverAD();
    else
        pSolver = BackslashSolverAD();
    end
    linsolve = CPRSolverAD('ellipticSolver', pSolver);
end

if opt.DynamicTimesteps
    timestepper = ...
    IterationCountTimeStepSelector('firstRampupStepRelative', 0.1, ...
                                   'firstRampupStep',         1*day);
end
solver = NonLinearSolver('timeStepSelector', timestepper, ...
                         'LinearSolver', linsolve, varg{:});
end
