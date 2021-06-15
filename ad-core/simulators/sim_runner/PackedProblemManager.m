classdef PackedProblemManager < handle
    % Class for managing a set of packed simulation problems. See
    % packSimulationProblem for more details.
    properties
        numThreadsPerSimulation = 1; % Number of threads in each simulation
        maxSimultaneousSimulations = maxNumCompThreads(); % Max number of simultaneous simulations in the background
        verbose = mrstVerbose(); % Verbosity
        displayProgressType % Type of progress update
        delay = 1; % Delay between polling (in seconds)
        showOnlyActive = []; % Show only active simulations in progress output
    end
    
    properties (Access = protected)
        packedProblems = {}; % Cell array of packed problems
        ncomplete = []; % Number of completed step
        nsteps = []; % Number of timesteps in schedule
        problem_is_executing = []; % Indicator if active
        N % Number of problems
        handle = struct('figure', [], 'text', [], 'iteration', 0); % Handle for monitorProgress
    end
    
    methods
        function ppm = PackedProblemManager(problems, varargin)
            ppm = merge_options(ppm, varargin{:});
            if isempty(ppm.displayProgressType)
                if usejava('desktop')
                    ppm.displayProgressType = 'graphical';
                else
                    ppm.displayProgressType = 'text';
                end
            end
            if ~iscell(problems)
                problems = {problems};
            end
            ppm.packedProblems = reshape(problems, [], 1);
            ppm.N = numel(problems);
            if isempty(ppm.showOnlyActive)
                ppm.showOnlyActive = ppm.N > 20;
            end
            ppm.nsteps = ppm.problemfun(@(x) numel(x.SimulatorSetup.schedule.step.val));
            ppm.problem_is_executing = false(ppm.N, 1);
        end
        % --- Simulator gateways
        function simulateProblemsBatch(ppm, index)
            % Simulate a number of problems in the background (index
            % optional)
            ppm.handle = struct('figure', [], 'text', [], 'iteration', 0);
            if nargin == 1
                index = ':';
            end
            done = false;
            while ~done
                [launched, sim_is_done] = ppm.launchBatchSimulations(index);
                done = all(sim_is_done);
                ppm.monitorProgress(true, index);
                if ppm.delay > 0
                    pause(ppm.delay);
                end
                if any(launched) && ppm.verbose > 0
                    fprintf('Launched %d new sessions...\n', numel(launched));
                end
            end
        end
        
        function simulateProblem(ppm, varargin)
            % Simulate one or more problems (directly in current session)
            if mod(nargin - 1, 2) == 1
                index = varargin{1};
                varargin = varargin(2:end);
            else
                index = ':';
            end
            simulatePackedProblem(ppm.packedProblems(index), varargin{:});
        end
        % --- Utilities
        function monitorProgress(ppm, single_update, index)
            % Show progress indicator for background simulation of problems
            if nargin < 3
                index = ':';
            end
            problems = ppm.packedProblems(index);
            if nargin == 1
                single_update = false;
            end
            indices = (1:ppm.N)';
            if ppm.showOnlyActive
                active = ppm.problem_is_executing(index);
                fh = ppm.handle.figure;
                if any(active) || (~isempty(fh) && ishandle(fh))
                    isDone = ppm.getNumberOfCompleteSteps(false) >= ppm.nsteps;
                    alreadyDone = sum(isDone(~active));
                else
                    % Display everything, for an overview
                    active = ':';
                    alreadyDone = 0;
                end
            else
                active = ':';
                alreadyDone = 0;
            end
            indices = indices(active);
            activeProblems = problems(active);
            switch lower(ppm.displayProgressType)
                case {'graphical', 'text', 'text-scrolling'}
                    useFigure = strcmpi(ppm.displayProgressType, 'graphical');
                    ppm.handle = monitorBackgroundSimulations(activeProblems, ...
                        'useFigure', useFigure, 'handle', ppm.handle, ...
                        'dynamicText', ~strcmpi(ppm.displayProgressType, 'text-scrolling'), ...
                        'totalNumberOfCases', numel(problems), 'totalProgress', alreadyDone, ...
                        'singleUpdate', single_update, 'indices', indices);
                otherwise

            end
        end
        
        function [launched, done] = launchBatchSimulations(ppm, index)
            % For a given index set, try to launch background simulations.
            % Will only launch simulators within capacity dictated by
            % maxSimultaneousSimulations property
            if nargin == 1 || ischar(index)
                index = (1:ppm.N)';
            end
            if islogical(index)
                index = find(index);
            end
            completed_steps = ppm.getNumberOfCompleteSteps(true);
            done = completed_steps(index) >= ppm.nsteps(index);
            % Check if some problems are already executing
            n = getNumberOfRunningSimulations(ppm);
            capacity = ppm.maxSimultaneousSimulations - n;
            if capacity > 0
                candidates = index(~done & ~ppm.problem_is_executing(index));
                next = candidates(1:min(capacity, numel(candidates)));
                for i = 1:numel(next)
                    ppm.launchBatchSimulationNoCapacityCheck(next(i));
                end
            else
                next = [];
            end
            launched = next;
        end
        
        function launchBatchSimulationNoCapacityCheck(ppm, index)
            % Launch a single simulation in the background (no check for
            % capacity! Use launchBatchSimulations for capacity check)
            if islogical(index)
                index = find(index);
            end
            assert(numel(index) == 1);
            assert(~ppm.problem_is_executing(index));
            nthread = ppm.numThreadsPerSimulation;
            if nthread == 1
                parg = '-singleCompThread';
            else
                parg = '';
            end
            simulateProblemBackground(ppm, index, 'extra_arg', parg, 'maxNumCompThreads', nthread);
            ppm.problem_is_executing(index) = true;
        end
        
        function n = getNumberOfRunningSimulations(ppm)
            % Get the number of current running simulations
            ppm.updateNumberOfRunningSimulations();
            n = sum(ppm.problem_is_executing);
        end
        
        function updateNumberOfRunningSimulations(ppm)
            % Update the currently running simulations (include checking if
            % the currently running simulations have completed)
            n = ppm.getNumberOfCompleteSteps(true);
            ppm.problem_is_executing = ppm.problem_is_executing & ~(n >= ppm.nsteps);
        end

        function simulateProblemBackground(ppm, varargin)
            % Simulate a problem in background (without interacting with
            % the capacity of the manager at all)
            if mod(nargin - 1, 2) == 1
                index = varargin{1};
                varargin = varargin(2:end);
            else
                index = (1:ppm.N)';
            end
            problems = ppm.packedProblems(index);
            for i = 1:numel(problems)
                simulatePackedProblemBackground(problems{i}, varargin{:});
            end
        end
        % --- Problem utilities
        function out = problemfun(ppm, fn, unif)
            if nargin < 3
                unif = true;
            end
            out = cellfun(fn, ppm.packedProblems, 'UniformOutput', unif);
        end
        
        function n = getNumberOfCompleteSteps(ppm, update)
            if nargin == 1
                update = true;
            end
            if update
                ppm.updateNumberOfCompletedSteps();
            end
            n = ppm.ncomplete;
        end

        function updateNumberOfCompletedSteps(ppm)
            ppm.ncomplete = ppm.problemfun(@(x) numelData(x.OutputHandlers.states));
        end
        
        function problems = getPackedProblems(ppm, index)
            problems = ppm.packedProblems;
            if nargin > 1
                problems = problems(index);
            end
        end
    end
end

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
