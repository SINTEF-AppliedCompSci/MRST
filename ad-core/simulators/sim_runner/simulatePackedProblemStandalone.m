function simulatePackedProblemStandalone(pth)
%Stand-alone solver for running packed problems programmatically
%
% SYNOPSIS:
%   simulatePackedProblemStandalone('/path/to/saved/problem')
%
% REQUIRED PARAMETERS:
%   pth - Path to saved struct in .mat format containing the fields
%         'modlist' (cell-array of all modules required to simulate) and
%         'problem' (saved simulation problem from packSimulationProblem).
%         Note that modules are stored seperately, to make it possible to
%         load modules before loading the problem itself, which requires
%         all classes to be on the path.
%
% RETURNS:
%   Nothing.
%
% SEE ALSO:
%   packSimulationProblem, simulatePackedProblemBackground

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

    p_path = fullfile(pth, 'problem.mat');
    % Load just the modules
    tmp = load(p_path, 'modlist');
    mrstModule('reset', tmp.modlist{:});
    % Now we can load the whole problem
    tmp = load(p_path, 'problem', 'opt');
    nthread = tmp.opt.maxNumCompThreads;
    if isfinite(nthread) && maxNumCompThreads() ~= nthread
        maxNumCompThreads(nthread);
    end
    lockpath = fullfile(pth, 'lock.mrst');
    if exist(lockpath, 'dir')
        error('This simulation is already running!')
    else
        fclose(fopen(lockpath, 'w'));
        try
            simulatePackedProblem(tmp.problem);
        catch

        end
        delete(lockpath);
    end
end
