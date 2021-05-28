function info = simulatePackedProblemBackground(problem, varargin)
%Simulate a packed simulation problem as a seperate Matlab thread
%
% SYNOPSIS:
%   info = simulatePackedProblemBackground(problem)
%
% REQUIRED PARAMETERS:
%   problem - Packed simulation problem from 'packSimulationProblem'.
%
% OPTIONAL PARAMETERS:
%   various - Parameters for setting up the seperate Matlab session.
%
% RETURNS:
%   info    - Struct containing information about where the log, lock and
%             saved problems are stored.
% NOTE:
%   This routine generates an additional Matlab session, either running in
%   the background (Linux / Mac OS) or as a seperate Window. This allows
%   for parallel running of simulations with just a base Matlab license,
%   but has some limitations / caveats:
%      - The spawned session will run until completion or error! If you
%      start a very large simulation and want to stop the simulation, you
%      must terminate the process.
%      - You must have Matlab on the path so that system('matlab') starts a
%      Matlab session. This is normally the case for standard
%      installations.
%
% SEE ALSO:
%   simulatePackedProblem, packSimulationProblem

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

    opt = struct('workdir', fullfile(mrstOutputDirectory(), 'sim_runner'), ...
                 'linux_arg', '-nodisplay', ...
                 'win_arg', '', ...
                 'matlab_arg', '-nosplash -noFigureWindows -nodesktop', ...
                 'extra_arg', '', ...
                 'verbose', mrstVerbose(), ...
                 'maxNumCompThreads', inf);
    opt = merge_options(opt, varargin{:});
    basepath = fullfile(opt.workdir, [problem.BaseName, '_', problem.Name]);
    dispif(opt.verbose, 'Storing problem %s: %s to %s...', problem.BaseName, problem.Name, basepath);
    if not(exist(basepath, 'dir'))
        mkdir(basepath);
    end
    fn = makepath(basepath);
    logfile = makepath(fn, 'simulation.log');
    modlist = problem.Modules;
    pathmap = capture_mrstpath();
    pathmap_file = makepath(fn, 'pathmap.mat');
    save(pathmap_file, 'pathmap');
    save(makepath(fn, 'problem'), 'problem', 'modlist', 'opt');
    matlab_options = build_matlab_options(opt);
    dispif(opt.verbose, 'Ok!\n');
    dispif(opt.verbose, 'Initializing background Matlab session...');
    str = sprintf('matlab %s -sd "%s" -logfile "%s" -r "pm = load(''%s''); mrstPath(''reregister'', pm.pathmap{:}); clear pm; mrstModule add ad-core; simulatePackedProblemStandalone(''%s''); exit()" &', ...
        matlab_options, ...% Matlab options
        makepath(ROOTDIR()), ...% Startup dir - MRST root
        logfile, ... % Log file path
        pathmap_file, ...
        fn ... % Problem path
        );
    system(str);
    dispif(opt.verbose, ' Ok! Simulation launched.\n');
    info = struct('path', basepath);
end

function options = build_matlab_options(opt)
    if ispc()
        % Windows arguments
        base = opt.win_arg;
    else
        % Mac OS and linux uses same parameters
        base = opt.linux_arg;
    end
    options = sprintf('%s %s %s', base, opt.matlab_arg, opt.extra_arg);
end

function str = makepath(varargin)
    str = fullfile(varargin{:});
    seperator = filesep();
    if strcmp(seperator, '\')
        str = strrep(str, '\', '\\');
    end
end

function pathmap = capture_mrstpath()
    mods  = mrstPath();
    paths = mrstPath(mods{:});

    pathmap = [ reshape(mods , 1, []) ; ...
                reshape(paths, 1, []) ];
end
