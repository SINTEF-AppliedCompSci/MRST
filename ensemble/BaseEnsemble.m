classdef BaseEnsemble < handle
    % Base class that facilitates ensembles in MRST.
    %
    % SYNOPSIS:
    %   ensembles = Ensemble(samples, 'name', name, ...)
    %
    % DESCRIPTION:
    %   This class is used to organize, set-up, run and post-process 
    %   ensemble simulations in MRST. It can be used for a varity of 
    %   problems, such as uncertainty quantification, optimization under 
    %   uncertainty, and history matching.
    % 
    % REQUIRED PARAMETERS:
    %
    %   samples - Sample object, typically of a superclass of BaseSamples.
    %             It represents stochastic parameters or configurations
    %             that are specific for each individual ensemble member.
    %             This object will also have functionality to map sample
    %             realizations to the base problem and thereby create the
    %             ensemble member simulations.
    %
    % OPTIONAL PARAMETERS:
    %   'name'      - Ensemble name. If not provided and 'setup' is provided, 
    %                 setup.name will be used.
    %   'setup'     - Common settings for all ensemble members, e.g., an
    %                 MRSTExample or a class/struct containing 
    %                 getPackedSimulationProblem (see method getBaseProblem 
    %                 for usage). If left empty, sample-property must provide 
    %                 entire setup of sample problems (see method getSampleProblem) 
    %   'directory' - Path to where we store results when simulating
    %                 the ensemble. Automatically generated if not provided
    %                 Default:
    %                 mrstOutputDirectory()/ensemble/ensemble.name
    %   'reset' - Boolean (default: false). Will delete any old results
    %             existing in the 'directory' so that any simulation
    %             results related to the new ensemble will have to be
    %             regenerated.
    %   'solve' - Function handle for simulating/running an ensemble
    %             member.
    %             Default: @simulatePackedProblem(problem)
    %   'simulationStrategy' - String to control how to run the ensemble.
    %                          Acceptable values:
    %                       'serial'     - No parallelization (default)
    %                       'parallel'   - Run ensemble in parallel using 
    %                           the Parallel Computing Toolbox,
    %                       'background' - Run ensemble members by spawning
    %                           new matlab sessions in the background.
    %   'maxWorkers' - Maximum number of parallel workers to use for
    %                  processing the ensemble members (default is system
    %                  dependent)
    %   'verbose' - Boolean (default: true). Print some extra info
    %   'verboseSimulation' - Boolean (default: false). Print the output of
    %                         each ensemble member simulation   
    %   'matlabBinary' - If 'simulationStrategy' is set to 'background',
    %                    it is possible to provide a specific location of
    %                    the matlab binary (typically used if you have
    %                    several matlab installations on your system).
    %   'backgroundEvalFn' - Function name for running ensemble members in
    %                        the background. 
    %                        Default: @simulateEnsembleMembersStandalone
    %
    %   If the 'mrstExample' parameter is a string, it is possible to add
    %   extra parameters which will then be passed on to the MRSTExample
    %   class.
    %
    % RETURNS:
    %   Class instance.
    %
    % SEE ALSO:
    %   `MRSTEnsemble`
    
    properties
        name
        directory
        setup
        num
        samples 
        % Function handle for solving a given problem
        solve = @(problem) simulatePackedProblem(problem); 
        
        simulationStrategy = 'serial';
        maxWorkers = maxNumCompThreads();
        cluster = [];
        spmdEnsemble; % Used for 'simulationStrategy' = 'spmd'

        
        verbose = true
        verboseSimulation = false
        
        simulationStatus  % handler to fetch simulation status
        
        backgroundEvalFn = @simulateEnsembleMembersStandalone
        matlabBinary = ''  

        hasBaseProblem
        
        NonLinearSolver = [];
    end
        
    
    methods
        %-----------------------------------------------------------------%
        function ensemble = BaseEnsemble(samples, varargin)
            [ensemble, extra] = merge_options(ensemble, varargin{:});
            ensemble.setupName();
            ensemble.setUpDirectory();
            % if there is a setup, this will provide the base-problem
            ensemble.hasBaseProblem = ~isempty(ensemble.setup);
            % other options
            opt = struct('reset',              false, ...
                         'prepareSimulation',  true,  ...
                         'statusFolder',       ''); 
            opt = merge_options(opt, extra{:});
            
            % Set samples
            ensemble.samples = samples;
            ensemble.num     = samples.num;

            % Setup simulation status handler
            ensemble.setupSimulationStatusHandler(opt.statusFolder);
            
            % Delete existing results if requested
            if opt.reset
                ensemble.reset('prompt',            false, ...
                               'prepareSimulation', false);
            end
            
            % Prepare ensemble
            if opt.prepareSimulation
                ensemble.prepareEnsembleSimulation();
            end
            
            ensemble.checkSetup();
        end
        
        %-----------------------------------------------------------------%
        function flag = reset(ensemble, varargin)
            % Deletes any old results so that we can start the ensemble
            % simulation fra scratch.      
            flag = true;
            opt = struct('prompt', true, 'prepareSimulation', true);
            opt = merge_options(opt, varargin{:});
            dataPath = ensemble.getDataPath();
            if ~exist(dataPath, 'dir'), return, end
            if opt.prompt
                prompt = sprintf(['Delete all data for %s? (sample '    , ...
                                  'data will not be deleted) y/n [n]: '], ...
                                              ensemble.name     );
                if ~strcmpi(input(prompt, 's'), 'y')
                    fprintf('Ok, will not remove files.\n');
                    flag = false;
                    return
                end
            end
            % Delete status data and corresponding status folder
            ensemble.simulationStatus.resetData();

            % Delete ensemble file
            if exist(fullfile(dataPath, 'ensemble.mat'), 'file')
                delete(fullfile(dataPath, 'ensemble.mat'));
            end
            % Delete base problem folder
            if exist(fullfile(dataPath, 'baseProblem'), 'dir')
                rmdir(fullfile(dataPath, 'baseProblem'), 's');
            end
            list = ls(dataPath);
            % Delete sample output directories
            samp = folderRegexp(list, '\d+\s', 'match');
            samp = cellfun(@str2num, samp);
            for s = samp
                fname = fullfile(dataPath, num2str(s));
                if exist(fname, 'dir')
                    rmdir(fname, 's');
                end
            end
            % Delete log files (background simulations only)
            logs = folderRegexp(list, 'log_\d+_\d+.mat', 'match');
            for l = logs
                delete(fullfile(dataPath, l{1}));
            end
            % Prepare ensemble
            if opt.prepareSimulation
                ensemble.prepareEnsembleSimulation();
            end
        end
        
        %-----------------------------------------------------------------%
        function problem = getBaseProblem(ensemble)
            if ensemble.hasBaseProblem
                % Returns the base problem in its pure form, without using any
                % of the stochastic samples.
                problem = ensemble.setup.getPackedSimulationProblem( ...
                               'Directory'      , ensemble.directory      , ...
                               'Name'           , 'baseProblem'           , ...
                               'NonLinearSolver', ensemble.NonLinearSolver);
                if ~isempty(ensemble.samples.processProblemFn)
                    problem = ensemble.samples.processProblemFn(problem);
                end
            else
                error('Ensemble has not been set up with a base-problem');
            end
        end
        
        %-----------------------------------------------------------------%
        function problem = getSampleProblem(ensemble, seed, varargin)
            if ensemble.hasBaseProblem
                % Setup using base-problem
                opt = struct('sample', []);
                opt = merge_options(opt, varargin{:});
                baseProblem = ensemble.getBaseProblem();
                problem     = ensemble.samples.getSampleProblem(baseProblem, seed, opt.sample);
            else
                % Setup without base-problem (path for result-handler must 
                % be provided
                problem = ensemble.samples.getSampleProblem(ensemble.directory, seed);
            end
        end
        
        %-----------------------------------------------------------------%
        function ensemble = setSolver(ensemble, solve)
            % Updates the function used for simulating the individual
            % ensemble members.
            ensemble.solve = solve;
        end
        
        %-----------------------------------------------------------------%
        function dataPath = getDataPath(ensemble)
            dataPath = ensemble.directory;
        end
              
        %-----------------------------------------------------------------%
        function defaultPath = getDefaultPath(ensemble)
            defaultPath = fullfile(mrstOutputDirectory(), ensemble.name);
        end
        
        %-----------------------------------------------------------------%
        function flag = hasSimulationOutput(ensemble, range)
            % check if a seed or range of seeds has stored output
            ids = ensemble.simulationStatus.getValidIds();
            n = ensemble.num;
            if isinf(n)
                n = max(range);
                if ~isempty(ids), n = max(n, max(ids)); end
            end
            flag = false(n, 1);
            flag(ids) = true;
            if nargin >= 2
                flag = flag(range);
            end
        end
        
        %-----------------------------------------------------------------%
        function flag = getSimulationStatus(ensemble, range)
            % check if a seed or range of seeds has simulation output and 
            % in addition the simulation status
            % -1: failed simulation
            %  0: not simulated
            %  1: successful simulation
            if nargin <= 1
                range = (1:ensemble.num)';
            end
            flag = double(ensemble.hasSimulationOutput(range));
            if any(flag)
                for k = reshape(find(flag), 1, [])
                    if ~ensemble.simulationStatus{range(k)}.success
                        flag(k) = -1;
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        function printSimulationStatus(ensemble)
            f = ensemble.getSimulationStatus;
            fprintf('Found output from %d out of %d simulations, %d failed.\n', ...
                nnz(f), ensemble.num, nnz(f<0));
        end
        %-----------------------------------------------------------------%
        function h = getSimulationHandlers(ensemble, prefix, range)
            % Get simulation result handlers with given prefix 
            if nargin <= 2
                range = (1:ensemble.num)';
            elseif islogical(range)
                range = find(range);
            end
            h = cell(1, numel(range));
            for k = 1:numel(h)
                h{k} = ResultHandler('dataPrefix',    prefix, ...
                                     'writeToDisk',   true, ...
                                     'dataDirectory', ensemble.directory, ...
                                     'dataFolder',    num2str(range(k)), ...
                                     'cleardir',      false);
            end
            % do a simple sanity check that prefix matches existing output from
            % simulations
            flag = ensemble.hasSimulationOutput(range);
            if any(flag)
                i1 = find(flag, 1, 'first');
                if isempty(h{i1}.getValidIds)
                    warning(...
                    'Unexpected empty data for handlers with prefix ''%s'' encountered.\n%s', ...
                    prefix, 'Check that prefix matches output from simulation.')
                end
            end

        end  
        
        %-----------------------------------------------------------------%
        function varargout = simulateEnsembleMember(ensemble, seed, varargin)
            % Run simulation for a specific ensemble member.
            %
            % SYNOPSIS:
            %   ensemble.simulateEnsembleMember(seed);
            %
            % PARAMETERS:
            %   seed - Integer specifying which ensemble member to run.
            %
            opt = struct('sample', []);
            opt = merge_options(opt, varargin{:});
            doSolve       = ~ensemble.hasSimulationOutput(seed);
            outputProblem = nargout > 0;
            
            if ~doSolve
                % Simulation is done - nothing to do here!
                if ensemble.verbose
                    fprintf(['Simulation output for seed %d found on ', ...
                             'disk (skipping)\n'], seed               );
                end
                if outputProblem
                    varargout{1} = ensemble.getSampleProblem(seed);
                    varargout{2} = ensemble.simulationStatus{seed};
                end
            else
                problem = ensemble.getSampleProblem(seed, 'sample', opt.sample);
                % Solve problem
                status = struct('success', true, 'message', []);
                try
                    ensemble.solve(problem);
                catch me
                    status.success = false;
                    status.message = me;
                end
                ensemble.simulationStatus{seed} = status;
                if outputProblem
                    varargout{1} = problem;
                    varargout{2} = status;
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function range = simulateEnsembleMembers(ensemble, varargin)
            % Run simulations for a set of ensemble members
            %
            % SYNOPSIS:
            %   ensemble.simulateEnsembleMembers();
            %   ensemble.simulateEnsembleMembers(20:30);
            %
            % OPTIONAL PARAMETERS:
            %   'range' - Specifies which ensemble members to run either as
            %             a vector or scalar.
            %           vector: The seeds (ids) for the ensemble members 
            %               that will be simulated
            %           scalar: Number of ensemble members to simulate. If
            %               results for some ensemble members already exist,
            %               an dditional 'range' number of members will be 
            %               run.
            
            opt = struct('range'    , inf, ...
                         'batchSize', inf);
            [opt, extraOpt] = merge_options(opt, varargin{:});
            % Validate optional input
            has_range = ~all(isinf(opt.range));
            has_size  = ~isinf(opt.batchSize);
            assert(~(has_range && has_size), 'Cannot specify both range and batchSize');
            if ~has_range && ~has_size
                % Neither range nor batchSize given, simulate full ensemble
                assert(~isinf(ensemble.num), ...
                    ['Ensemble size not defined, please use '                     , ...
                     'ensemble.simulateEnsembleMembers(''range'', range) instead']);
                opt.range     = 1:ensemble.num;
                opt.batchSize = numel(opt.range);
            end
            if has_size
                % Translate batch size to simulation range
                ids = ensemble.qoi.ResultHandler.getValidIds();
                if isempty(ids), ids = 0; end
                if max(ids) + opt.batchSize > ensemble.num
                    warning(['Requested ensemble members plus already ' , ...
                             'computed ensemble members are more than ' , ...
                             'ensemble size. Proceeding by trying to '  , ...
                             'simulate the last %d ensemble members.'  ], ...
                             opt.range)
                    opt.batchSize = ensemble.num - max(ids);
                end
                opt.range = (1:opt.batchSize) + max(ids);
            end
            % Validate range
            assert(max(opt.range) <= ensemble.num, ...
                'Requested range is outside range of available ensemble members');
            % Distribute ensemble members among workers/background sessions
            n        = ceil(numel(opt.range)/ensemble.maxWorkers);
            rangePos = repmat(n, ensemble.maxWorkers, 1);
            extra    = sum(rangePos) - numel(opt.range);
            rangePos(end-extra+1:end) = rangePos(end-extra+1:end) - 1;
            rangePos = cumsum([0; rangePos]) + 1;
            % Simulate with appropriate strategy
            switch ensemble.simulationStrategy
                case 'serial'
                    ensemble.simulateEnsembleMembersSerial(opt.range, extraOpt{:});
                case 'parallel'
                    ensemble.simulateEnsembleMembersParallel(opt.range, rangePos, extraOpt{:});
                case 'background'
                    ensemble.simulateEnsembleMembersBackground(opt.range, rangePos, extraOpt{:});
                case 'spmd'
                    ensemble.simulateEnsembleMembersSPMD(opt.range, rangePos, extraOpt{:});
                otherwise
                    error(['Invalid simulationStrategy "', ensemble.simulationStrategy, '"'])
            end
            range = opt.range;
        end
        
        %-----------------------------------------------------------------%
        function simulateEnsembleMembersCore(ensemble, range)
            % Runs simulations according to the given range.
           range = reshape(range, 1, []);
           for i = 1:numel(range)
               seed = range(i);
               % Print info
               if ensemble.verbose
                   progress = floor(100*(i-1)/numel(range));
                   fprintf(['(%d%%)\tSimulating ensemble member %d ', ...
                            'among ensemble members %d to %d...\n' ], ...
                           progress, seed, range(1), range(end)    );
               end
               % Run simulation
               if ensemble.verboseSimulation
                   ensemble.simulateEnsembleMember(seed);
               else
                   evalc('ensemble.simulateEnsembleMember(seed)');
               end
           end      
           % Print info
           if ensemble.verbose
               fprintf(['(100%%)\tDone simulating ensemble members ', ...
                        '%d to %d \n'], range(1), range(end)        );
           end
        end
        
        %-----------------------------------------------------------------%
        function checkSetup(ensemble)
            s = ensemble.setup;
            if ~isempty(s)
                if isstruct(s)
                    ok = isfield(s, 'getPackedSimulationProblem');
                elseif isobject(s)
                    ok = ismethod(s, 'getPackedSimulationProblem');
                else 
                    error('Unexpected format of setup: %s', class(s));
                end
                if ~ok
                    error('Property ''status'' does not contain required %s', ...
                          'field/prop ''getPackedSimulationProblem''');
                end
            end
        end
        
    end % methods
    
    methods (Access = protected)
        
        %-----------------------------------------------------------------%
        function setupName(ensemble)
            if isempty(ensemble.name)
                % fetch name from setup 
                if ~isempty(ensemble.setup)
                    ensemble.name = ensemble.setup.name;
                else
                    error('Property ''name'' cannot be defaulted for BaseEnsemble');
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function setUpDirectory(ensemble)
            % if ensemble has setup, fetch
            
            % Set up directory name correctly.
            % This function is mainly provided so that sub-classes can
            % define other directory structures.
            if isempty(ensemble.directory)
                ensemble.directory = ensemble.getDefaultPath();
            end
        end
            
        %-----------------------------------------------------------------%
        function setupSimulationStatusHandler(ensemble, subFolder)
            if isempty(subFolder)
                % setup in main folder
                [ddir, dfolder] = fileparts(ensemble.directory);
            else
                [ddir, dfolder] = deal(ensemble.directory, subFolder);
            end
            ensemble.simulationStatus = ...
                ResultHandler('dataPrefix',    'simulationStatus', ...
                              'writeToDisk',   true, ...
                              'dataDirectory', ddir, ...
                              'dataFolder',    dfolder, ...
                              'cleardir',      false);
        end         
        %-----------------------------------------------------------------% 
        function prepareEnsembleSimulation(ensemble, varargin)
            % INTERNAL
            % Set parameters related to possible parallel execution scheme
            % for simulating the ensemble members. 
            % Called by the constructor.
            
            opt = struct('force', false);
            opt = merge_options(opt, varargin{:});
            
            if strcmpi(ensemble.simulationStrategy, 'serial')
                warning(['Serial ensemble simulations will take a '     , ...
                         'long time, and should only be used for '      , ...
                         'debugging or small ensembles. Consider '      , ...
                         'using ''background'', ''parallel'' or ''spmd'' instead.'])
            end
            
            if any(strcmpi(ensemble.simulationStrategy, {'parallel', 'spmd'})) 
                % Check if parallel toolbox is availiable
                if isempty(ver('parallel'))
                    warning(['Parallel computing toolbox not available. ', ...
                             'Switching to background simulation']        );
                    ensemble.simulationStrategy = 'background';
                end
            end
            
            if any(strcmpi(ensemble.simulationStrategy, {'parallel', 'background'}))
                fn = fullfile(ensemble.getDataPath(), 'ensemble.mat');
                if ~exist(fn, 'file') || opt.force
                    if isprop(ensemble, 'figures')
                        for f = fieldnames(ensemble.figures)'
                            if isgraphics(ensemble.figures.(f{1}))
                                delete(ensemble.figures.(f{1}));
                            end
                            ensemble.figures.(f{1}) = [];
                        end
                    end
                    save(fn, 'ensemble');
                end
            end
            
            if strcmpi(ensemble.simulationStrategy, {'parallel'})
                if isempty(ensemble.cluster)
                    ensemble.cluster = parcluster('local');
                end
                ensemble.cluster.Jobs.delete();
            end
            
            if strcmpi(ensemble.simulationStrategy, 'spmd')
                 % Use parallel toolbox. Check if we have started a
                % parallel session already, and whether it has the
                % correct number of workers
                if ensemble.maxWorkers > maxNumCompThreads()
                    warning(['Requested number of workes is greater ' , ...
                             'than maxNumCompThreads (%d) reducing.'  ], ...
                             maxNumCompThreads()                       );
                    ensemble.maxWorkers = maxNumCompThreads();
                end
                if isempty(gcp('nocreate'))
                    parpool(ensemble.maxWorkers);
                elseif gcp('nocreate').NumWorkers ~= ensemble.maxWorkers
                    delete(gcp);
                    parpool(ensemble.maxWorkers);
                end

                if ~all(exist(ensemble.spmdEnsemble)) || opt.force %#ok
                    % Communicate ensemble to all workers
                    spmd
                        spmdEns = ensemble;
                    end
                    ensemble.spmdEnsemble = spmdEns;
                end
                
            end
        end
        
        %-----------------------------------------------------------------%
        function simulateEnsembleMembersSerial(ensemble, range, varargin)
            % Runs simulation for the given range in a serial configuration
            % (no parallel processing).
            ensemble.simulateEnsembleMembersCore(range);
        end
        
        %-----------------------------------------------------------------%
        function simulateEnsembleMembersParallel(ensemble, range, rangePos, varargin)
            % Runs simulation of ensemble members according to range across
            % parallel workers. The subset of range given to each worker is
            % specified by rangePos.
            %
            % NOTE:
            %   This function requires the Parallel Computing Toolbox.
            
            opt = struct('plotProgress', true, ...
                         'plotIntermediateQoI', true);
            [opt, extraOpt] = merge_options(opt, varargin{:});
            
            % Check that the parpool and spmdEnsemble is valid
            ensemble.prepareEnsembleSimulation()
            fileName = fullfile(ensemble.getDataPath(), 'ensemble.mat');
            n = 0;
            job = cell(ensemble.maxWorkers(), 1);
            for i = 1:numel(rangePos)-1
                r = range(rangePos(i):rangePos(i+1)-1);
                if isempty(r), continue; end
                n = n+1;
                job{i} = batch(ensemble.cluster, ensemble.backgroundEvalFn, 0, {fileName, r}, ...
                                                    'AttachedFiles', fileName);
            end
            fprintf(['Started %d new Matlab batch jobs. ', ...
                     'Waiting for simulations ...\n'    ], n);
            if ismethod(ensemble, 'plotProgress')
                %[varargout{1}, varargout{2}] = ensemble.plotProgress(range);
                ensemble.plotProgress(range, ...
                                      opt.plotProgress, opt.plotIntermediateQoI, ...
                                      extraOpt{:});
            end
            cellfun(@(job) delete(job), job); clear job;
        end
        
        %-----------------------------------------------------------------%
        function simulateEnsembleMembersBackground(ensemble, range, rangePos, varargin)
            % Runs simulation of ensemble members according to range by
            % launching matlab sessions in the background. Each matlab
            % session is responsible of running a subset of the ensemble
            % members specified in range according to rangePos.
            
            opt = struct('plotProgress', true, ...
                         'plotIntermediateQoI', true);
            [opt, extraOpt] = merge_options(opt, varargin{:});
            % Get path to mat file hodling the ensemble
            fileName = fullfile(ensemble.getDataPath(), 'ensemble.mat');
            % Make progress filename
            progressFileNm = @(r) fullfile(ensemble.getDataPath(), ...
                        ['log', '_', num2str(r(1)), '_', num2str(r(end)), '.mat']);
            % Get active MRST modules
            moduleList = mrstModule();
            n = 0; % Counter number of spawned sessions
            for i = 1:numel(rangePos)-1
                % Get local range
                r = range(rangePos(i):rangePos(i+1)-1);
                if isempty(r), continue; end
                n = n+1;
                % Spawn new MATLAB session
                evalFunWrapper(ensemble.backgroundEvalFn, {fileName, r}, ...
                               'progressFileNm', progressFileNm(r)     , ...
                               'moduleList'    , moduleList            , ...
                               'matlabBinary'  , ensemble.matlabBinary );
            end
            fprintf(['Started %d new Matlab sessions. ', ...
                     'Waiting for simulations ...\n'  ], n);
            if ismethod(ensemble, 'plotProgress')
                ensemble.plotProgress(range, ...
                                      opt.plotProgress, opt.plotIntermediateQoI, ...
                                      extraOpt{:});
            end
        end  
        
        %-----------------------------------------------------------------%
        function simulateEnsembleMembersSPMD(ensemble, range, rangePos, varargin)
            % Runs simulation of ensemble members according to range across
            % parallel workers. The subset of range given to each worker is
            % specified by rangePos.
            %
            % NOTE:
            %   This function requires the Parallel Computing Toolbox.
            
            % Check that the parpool and spmdEnsemble is valid
            ensemble.prepareEnsembleSimulation()
            
            spmdEns = ensemble.spmdEnsemble;
            spmd
                spmdRange = range(rangePos(labindex):rangePos(labindex+1)-1);
                spmdEns.simulateEnsembleMembersCore(spmdRange);
            end
            fprintf('simulateEnsembleMembersParallel completed\n');
        end
    end
end

%% Helpers
function matches = folderRegexp(list, expression, outputFormat)
    if size(list, 1) > 1
        % Windows behavior
        % Pad with a space at the end and reshape into a single
        % long string.
        pad = repmat(' ', size(list, 1), 1);
        list = reshape([list, pad]', 1, []);
    end
    matches = regexp(list, expression, outputFormat);
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
