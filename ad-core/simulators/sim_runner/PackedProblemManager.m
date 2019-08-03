classdef PackedProblemManager < handle
    % Class for managing a set of packed simulation problems. See
    % packSimulationProblem for more details.
    properties
        numThreadsPerSimulation = 1; % Number of threads in each simulation
        maxSimultaneousSimulations = maxNumCompThreads(); % Max number of simultaneous simulations in the background
        verbose = mrstVerbose(); % Verbosity
        displayProgressType % Type of progress update
        delay = 1; % Delay between polling (in seconds)
    end
    
    properties (Access = protected)
        packedProblems = {};
        ncomplete = [];
        nsteps = [];
        problem_is_executing = [];
        N
        handle = struct('figure', [], 'text', []);
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
            ppm.nsteps = ppm.problemfun(@(x) numel(x.SimulatorSetup.schedule.step.val));
            ppm.problem_is_executing = false(ppm.N, 1);
        end
        % --- Simulator gateways
        function simulateProblemsBatch(ppm, index)
            ppm.handle = struct('figure', [], 'text', []);
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
                    fprintf('Launched %d new sessions...\n', sum(launched));
                end
            end
        end
        
        function monitorProgress(ppm, single_update, index)
            problems = ppm.packedProblems;
            if nargin > 2
                problems = ppm.packedProblems(index);
            end
            if nargin == 1
                single_update = false;
            end
            
            switch lower(ppm.displayProgressType)
                case {'graphical', 'text'}
                    useFigure = strcmpi(ppm.displayProgressType, 'graphical');
                    ppm.handle = monitorBackgroundSimulations(problems, ...
                        'useFigure', useFigure, 'handle', ppm.handle, 'singleUpdate', single_update);
                otherwise

            end
        end
        
        function [launched, done] = launchBatchSimulations(ppm, index)
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
                    ppm.launchBatchSimulation(next(i));
                end
            else
                next = [];
            end
            launched = next;
        end
        
        function launchBatchSimulation(ppm, index)
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
            ppm.updateNumberOfRunningSimulations();
            n = sum(ppm.problem_is_executing);
        end
        
        function updateNumberOfRunningSimulations(ppm)
            n = ppm.getNumberOfCompleteSteps(true);
            ppm.problem_is_executing = ppm.problem_is_executing & ~(n >= ppm.nsteps);
        end
        
        function simulateProblem(ppm, varargin)
            if mod(nargin - 1, 2) == 1
                index = varargin{1};
                varargin = varargin(2:end);
            else
                index = ':';
            end
            simulatePackedProblem(ppm.packedProblems(index), varargin{:});
        end

        function simulateProblemBackground(ppm, varargin)
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
    end
end