function [simRes] = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, varargin)
% runSchedule -- Run simulation based on schedule.
%
% SYNOPSIS:
%   simRes = runSchedule(resSolInit, G, S, W, fluid, schedule, pn, pv, ...)
%
% DESCRIPTION:
%
% PARAMETERS:
%   resSolInit  -
%   G           - Grid data structure.
%   S           -
%   W           -
%   fluid       -
%   schedule    -
%
%
% RETURNS:
%   simRes      - (numSteps+1) x 1 structure having fields
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
                 'VerboseLevel', 2);
opt     = merge_options(opt, varargin{:});

verboseLevel2 = opt.Verbose || (opt.VerboseLevel == 2);
verboseLevel1 = opt.Verbose || (opt.VerboseLevel > 0);

numSteps = numel(schedule);
resSol  = resSolInit;

% Initial conditions
simRes(1).timeInterval  = [0 0];
simRes(1).resSol        = resSol;
simRes(1).wellSol       = resSol.wellSol;
if verboseLevel2, dispSchedule(schedule); end

% Pick solver
if strcmp(S.type, 'hybrid')
   solver = 'hybrid';
else
   solver = 'mixed';
end

if verboseLevel1, fprintf('\n******* Starting forward simulation *******\n'); end
for k = 1 : numSteps
    if verboseLevel1, fprintf('Time step %3d of %3d,   ', k, numSteps); end
    W     = updateWells(W, schedule(k));
    interval = schedule(k).timeInterval;
    dt       = interval(2) - interval(1);

    % ---- Pressure Equation -----
    if verboseLevel1, fprintf('Pressure:'); tic; end
    resSol = solveIncompFlowLocal(resSol, G, S, fluid, ...
                                        'wells', W, 'Solver', solver);

    if verboseLevel1, t = toc; fprintf('%9.3f sec,   ', t); end

    % ---- Saturation Equation ---
    if verboseLevel1, fprintf('Transport:'); tic; end

    resSol = implicitTransport(resSol, G, dt, rock, fluid,  ...
                               'wells', W,                           ...
                               'nltol', 1.0e-6, 'lstrials', 50,      ...
                               'maxnewt', 100,  'tsref',  15,        ...
                               'verbose', opt.Verbose);

    if verboseLevel1, t = toc; fprintf('%9.3f sec\n', t); end

    % update simRes structure
    simRes(k+1).timeInterval  = interval;
    simRes(k+1).resSol        = resSol;
    simRes(k+1).wellSol       = resSol.wellSol;
end


