function [adjRes] = runAdjoint(simRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, varargin)
% runAdjoint -- Run adjoint simulation based on simRes and schedule.
%
% SYNOPSIS:
%   adjRes = runAdjoint(simRes, G, S, W, fluid, schedule, objective, pn, pv, ...)
%
% DESCRIPTION:
%
% PARAMETERS:
%   simRes      -
%   G           - Grid data structure.
%   S           -
%   W           -
%   fluid       -
%   schedule    -
%   objective   - function handle
%
%
% RETURNS:
%   adjRes      - numSteps x 1 structure having fields
%                   - timeInterval
%                   - resSol
%                   - wellSol
%
%
% SEE ALSO:

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


opt     = struct('Verbose',  mrstVerbose , ...
                 'VerboseLevel', 0);
opt     = merge_options(opt, varargin{:});
verboseLevel2 = opt.Verbose || (opt.VerboseLevel == 2);
verboseLevel1 = opt.Verbose || (opt.VerboseLevel > 0);

numSteps = numel(schedule);
adjRes   = [];
obj      = objectiveFunction(G, S, W, rock, fluid, simRes, schedule, controls);
if verboseLevel1, fprintf('\n******* Starting adjoint simulation *******\n'); end
for k = numSteps : -1 : 1
    if verboseLevel1, fprintf('Time step %3d of %3d,   ', k, numSteps); end
    W    = updateWells(W, schedule(k));

     % ---- Analogue of Saturation Equation
    if verboseLevel1, fprintf('Transport:'); tic; end
    adjRes = solveAdjointTransportSystem(G, S, W, rock, fluid, ...
                                         simRes, adjRes, obj);
    if verboseLevel1, t = toc; fprintf('%9.3f sec,   ', t); end

    % ---- Analogue of Pressure Equation
    if verboseLevel1, fprintf('Pressure:'); tic; end
    adjRes = solveAdjointPressureSystem(G, S, W, rock, fluid,  ...
                                        simRes, adjRes, obj);
    if verboseLevel1, t = toc; fprintf('%9.3f sec\n', t); end
end
