classdef MRSTEnsemble
    % Class that facilitates ensembles in MRST.
    %
    % SYNOPSIS:
    %   ensembles = Ensemble(mrstExample, samples, qoi)
    %
    % DESCRIPTION:
    %   This class is used to organize, set-up, run and post-process 
    %   ensemble simulations in MRST. It can be used for a varity of 
    %   problems, such as uncertainty quantification, optimization, and
    %   history matching.
    % 
    % REQUIRED PARAMETERS:
    %   mrstExample - Instance of MRSTExample, or a function name that 
    %                 generates a MRSTExample, that serves as a base 
    %                 problem for the ensemble. If will contain all aspects
    %                 and configurations that will be common for all
    %                 ensemble members
    %
    %   samples - Sample object, typically of a superclass of BaseSamples.
    %             It represents stochastic parameters or configurations
    %             that are specific for each individual ensemble member.
    %             This object will also have functionality to map sample
    %             realizations to the base problem and thereby create the
    %             ensemble member simulations.
    %
    %   qoi - Quantity of interest object, typically a superclass of
    %         BaseQoI. It contains a post-processing mapping of a simulated
    %         problem (ensemble member) to a quantity of interest, and is
    %         therefore used to store and read relevant simulation results
    %         for the ensemble.
    %
    % OPTIONAL PARAMETERS:
    %   'directory' - Path to where we store results when simulating
    %                 the ensemble. Automatically generated if not provided
    %                 Default:
    %                 mrstOutputDirectory()/ensemble/mrstExample.name
    %   'reset' - Boolean (default: false). Will delete any old results
    %             existing in the 'directory' so that any simulation
    %             results related to the new ensemble will have to be
    %             regenerated.
    %   'storeOutput' - Boolean. If false, the simulation results (states,
    %                   wellSols, etc) of individual ensemble members will 
    %                   be deleted after the quantity of interest is 
    %                   computed and stored.
    %                   Default: false
    %   'solve' - Function handle for simulating/running an ensemble
    %             member.
    %             Default: @simulatePackedProblem(problem)
    %   'simulationType' - String to control how to run the ensemble.
    %                      Acceptable values:
    %                      'serial'     - No parallelization (default)
    %                      'parallel'   - Run ensemble in parallel using 
    %                          the Parallel Computing Toolbox,
    %                      'background' - Run ensemble members by spawning
    %                         new matlab sessions in the background.
    %   'maxWorkers' - Maximum number of parallel workers to use for
    %                  processing the ensemble members (default is system
    %                  dependent)
    %   'verbose' - Boolean (default: true). Print some extra info
    %   'verboseSimulation' - Boolean (default: false). Print the output of
    %                         each ensemble member simulation   
    %   'matlabBinary' - If the 'simulationType' is set to 'background',
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
    %   `MRSTExample`, `BaseSamples`, `BaseQoI`
    
    properties
        setup
        num
        samples
        qoi
        
        directory 
        storeOutput = false
        
        
        % Function handle for solving a given problem
        solve = @(problem) simulatePackedProblem(problem); 
        
        simulationType = 'serial';
        maxWorkers = maxNumCompThreads();
        spmdEnsemble; % Used for 'simulationType' = 'parallel'
        
        verbose = true
        verboseSimulation = false
        
        backgroundEvalFn = @simulateEnsembleMembersStandalone
        matlabBinary = ''
        
    end
    
    methods
        %-----------------------------------------------------------------%
        function ensemble = MRSTEnsemble(mrstExample, samples, qoi, varargin)
            [ensemble, extra] = merge_options(ensemble, varargin{:});
            opt          = struct('reset', false);
            [opt, extra] = merge_options(opt, extra{:});
            % Set example. This defines the base problem
            if isa(mrstExample, 'MRSTExample')
                % Example given
                ensemble.setup = mrstExample;
            else
                % Example name given - set up example
                ensemble.setup = MRSTExample(mrstExample, extra{:});
            end
            % Set up directory
            if isempty(ensemble.directory)
                ensemble.directory = fullfile(mrstOutputDirectory(), ...
                                              'ensemble', ensemble.setup.name);
            end
            % Set samples
            ensemble.samples = samples;
            ensemble.num     = samples.num;
            % Validate qoi
            baseProblem  = ensemble.getBaseProblem();
            ensemble.qoi = qoi.validateQoI(baseProblem);
            % Delete existing results if requested
            if opt.reset
                ensemble.reset('prompt', false);
            end
            % Prepare ensemble
            ensemble = ensemble.prepareEnsemble();
        end
        
        %-----------------------------------------------------------------%
        function ensemble = prepareEnsemble(ensemble)
            % Set parameters related to possible parallel execution scheme
            % for simulating the ensemble members. 
            % Called by the constructor.
            
            if strcmpi(ensemble.simulationType, 'serial')
                warning(['Serial ensemble simulations will take a '     , ...
                         'long time, and should only be used for '      , ...
                         'debugging or small ensembles. Consider '      , ...
                         'using ''background'' or ''parallel'' instead.'])
            end
            
            if strcmpi(ensemble.simulationType, 'parallel')
                % Check if parallel toolbox is avialiable
                if isempty(ver('parallel'))
                    warning(['Parallel computing toolbox not aviailable. ', ...
                             'Switching to background simulation']        );
                    ensemble.simulationType = 'background';
                end
            end
            
            switch ensemble.simulationType
                case 'parallel'
                    % Use parallel toolbox. Check if we have started a
                    % parallel session already, and whether it has the
                    % correct number of workers
                    if isempty(gcp('nocreate'))
                        parpool(ensemble.maxWorkers);
                    elseif gcp('nocreate').NumWorkers ~= ensemble.maxWorkers
                        delete(gcp);
                        parpool(ensemble.maxWorkers);
                    end
                    % Communicate ensemble to all workers
                    spmd
                        spmdEns = ensemble;
                    end
                    ensemble.spmdEnsemble = spmdEns;
                case 'background'
                    % Run simulations in background sessions
                    fn = fullfile(ensemble.getDataPath(), 'ensemble.mat');
                    if ~exist(fn, 'file')
                        save(fn, 'ensemble');
                    end
            end
        end
        
        %-----------------------------------------------------------------%
        function ensemble = reset(ensemble, varargin)
            % Deletes any old results so that we can start the ensemble
            % simulation fra scratch.
            
            opt = struct('prompt', true);
            opt = merge_options(opt, varargin{:});
            dataPath = ensemble.getDataPath();
            if ~exist(dataPath, 'dir'), return, end
            if opt.prompt
                prompt = sprintf(['Delete all data for %s? (sample '    , ...
                                  'data will not be deleted) y/n [n]: '], ...
                                              ensemble.setup.name     );
                if ~strcmpi(input(prompt, 's'), 'y')
                    fprintf('Ok, will not remove files.\n');
                    return
                end
            end
            % Delete QoIs
            ensemble.qoi.ResultHandler.resetData();
            % Delet base problem folder
            if exist(fullfile(dataPath, 'baseProblem'), 'dir')
                rmdir(fullfile(dataPath, 'baseProblem'), 's');
            end
            % Delete ensemble file
            if exist(fullfile(dataPath, 'ensemble.mat'), 'file')
                delete(fullfile(dataPath, 'ensemble.mat'));
            end
            list = ls(dataPath);
            % Delete sample output directories
            samp = folderRegexp(list, '\d+\s', 'match');
            samp = cellfun(@str2num, samp);
            for s = samp
                rmdir(fullfile(dataPath, num2str(s)), 's');
            end
            % Delete log files (background simulations only)
            logs = folderRegexp(list, 'log_\d+_\d+.mat', 'match');
            for l = logs
                delete(fullfile(dataPath, l{1}));
            end
            % Prepare ensemble
            ensemble = ensemble.prepareEnsemble();
        end
        
        %-----------------------------------------------------------------%
        function problem = getBaseProblem(ensemble)
            % Returns the base problem in its pure form, without using any
            % of the stochastic samples.
            
            problem = ensemble.setup.getPackedSimulationProblem('Directory', ensemble.directory, 'Name', 'baseProblem');
        end
        
        %-----------------------------------------------------------------%
        function ensemble = setSolver(ensemble, solve)
            % Updates the function used for simulating the individual
            % ensemble members.
            
            ensemble.solve = solve;
        end
        
        %-----------------------------------------------------------------%
        function dataPath = getDataPath(ensemble)
            dataPath = ensemble.directory();
        end
        
        %-----------------------------------------------------------------%
        function simulateEnsembleMember(ensemble, seed, varargin)
            % Run simulation for a specific ensemble member.
            %
            % SYNOPSIS:
            %   ensemble.simulateEnsembleMember(seed);
            %
            % PARAMETERS:
            %   seed - Integer specifying which ensemble member to run.
            %
            
            if ensemble.qoi.isComputed(seed)
                % QoI is already computed - nothing to do here!
                if ensemble.verbose
                    disp(strcat("Simulation for ", num2str(seed), " found on disk"));
                end
                return;
            end
            % Get base problem
            baseProblem = ensemble.getBaseProblem;
            % Set up sample problem from seed
            problem = ensemble.samples.getSampleProblem(baseProblem, seed);
            % Solve problem
            ensemble.solve(problem);
            % Compute QoI
            ensemble.qoi.getQoI(problem);
            % Clear output if requested
            if ~ensemble.storeOutput
                clearPackedSimulatorOutput(problem, 'prompt', false);
            end
        end
        
        %-----------------------------------------------------------------%
        function simulateEnsembleMembers(ensemble, range, varargin)
            % Run simulations for a set of ensemble members
            %
            % SYNOPSIS:
            %   ensemble.simulateEnsembleMembers(range);
            %
            % PARAMETERS:
            %   range - Specifies which ensemble members to run either as a
            %           vector or scalar.
            %       vector: The seeds (ids) for the ensemble members that
            %           will be simulated
            %       scalar: Number of ensemble members to simulate. If
            %           results for some ensemble members already exist, an
            %           additional 'range' number of members will be run.
            if isscalar(range)
                ids = ensemble.qoi.ResultHandler.getValidIds();
                if isempty(ids), ids = 0; end
                range = (1:range) + max(ids);
            end
            n        = ceil(numel(range)/ensemble.maxWorkers);
            rangePos = repmat(n, ensemble.maxWorkers, 1);
            extra    = sum(rangePos) - numel(range);
            rangePos(end-extra+1:end) = rangePos(end-extra+1:end) - 1;
            rangePos = cumsum([0; rangePos]) + 1;
            switch ensemble.simulationType
                case 'serial'
                    ensemble.simulateEnsembleMembersSerial(range, rangePos, varargin{:});
                case 'parallel'
                    ensemble.simulateEnsembleMembersParallel(range, rangePos);
                case 'background'
                    ensemble.simulateEnsembleMembersBackground(range, rangePos, varargin{:});
            end
        end
        
        %-----------------------------------------------------------------%
        function simulateAllEnsembleMembers(ensemble, varargin)
            % Run simulations for all ensemble members
            %
            % SYNOPSIS:
            %   ensemble.simulateAllEnsembleMembers();
            %
            % NOTE:
            %   Does not work if stochastic samples are generated on the
            %   fly with a generatorFn. If so, see
            %   `ensemble.simulateEnsembleMembers(range)`
            assert(~isinf(ensemble.num), ...
                'Ensemble size not define, please use ensemble.simulateEnsembleMembers(range) instead');
            ensemble.simulateEnsembleMembers(1:ensemble.num, varargin{:});
        end
        
        %-----------------------------------------------------------------%
        function simulateEnsembleMembersCore(ensemble, range)
            % INTERNAL
            % Runs simulations according to the given range.
            
           range = reshape(range, 1, []);
           for i = 1:numel(range)
               seed = range(i);
               
               % Print info
               if ensemble.verbose
                   progress = floor(100*(i-1)/numel(range));
                   fprintf('(%d%%)\tSimulating ensemble member %d among ensemble members %d to %d...\n', ...
                           progress, seed, range(1), range(end))                    
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
               fprintf('(100%%)\tDone simulating ensemble members %d to %d \n', range(1), range(end));
           end
        end

        %-----------------------------------------------------------------%
        function simulateEnsembleMembersSerial(ensemble, range, rangePos, varargin)
            % INTERNAL 
            % Runs simulation for the given range in a serial configuration
            % (no parallel processing).
            ensemble.simulateEnsembleMembersCore(range);
        end
        
        %-----------------------------------------------------------------%
        function simulateEnsembleMembersParallel(ensemble, range, rangePos)
            % INTERNAL
            % Runs simulation of ensemble members according to range across
            % parallel workers. The subset of range given to each worker is
            % specified by rangePos.
            %
            % NOTE:
            %   This function requires the Parallel Computing Toolbox.
            spmdEns = ensemble.spmdEnsemble;
            spmd
                spmdRange = range(rangePos(labindex):rangePos(labindex+1)-1);
                spmdEns.simulateEnsembleMembersCore(spmdRange);
            end
            fprintf('simulateEnsembleMembersParallel completed\n');
        end
        
        %-----------------------------------------------------------------%
        function simulateEnsembleMembersBackground(ensemble, range, rangePos, varargin)
            % INTERNAL
            % Runs simulation of ensemble members according to range by
            % launching matlab sessions in the background. Each matlab
            % session is responsible of running a subset of the ensemble
            % members specified in range according to rangePos.
            
            opt = struct('plotProgress', false);
            opt = merge_options(opt, varargin{:});
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
                evalFunWrapper(ensemble.backgroundEvalFn, {fileName, r} , ...
                               'progressFileNm', progressFileNm(r)      , ...
                               'moduleList'    , moduleList             , ...
                               'matlabBinary'  , ensemble.matlabBinary);
            end
            fprintf(['Started %d new Matlab sessions. ', ...
                     'Waiting for simulations ...\n'  ], n);
            if opt.plotProgress
                ensemble.plotProgress(range);
            end
        end
        
        %-----------------------------------------------------------------%
        function progress = getEnsembleMemberProgress(ensemble, range)
            % INTERNAL 
            % Utility function for monitoring the progression of an
            % ensemble member that is being run right now.
            if nargin < 2, range = ensemble.num; end
            progress = zeros(numel(range),1);
            nsteps   = numel(ensemble.setup.schedule.step.val);
            for i = 1:numel(range)
                if exist(fullfile(ensemble.directory(), ...
                        [ensemble.qoi.ResultHandler.dataPrefix, num2str(range(i)), '.mat']), 'file')
                    progress(i) = inf;
                    continue
                end
                dataDir = fullfile(ensemble.directory(), num2str(range(i)));
                if ~exist(dataDir, 'dir')
                    continue;
                end
                files = ls(dataDir);
                progress(i) = numel(folderRegexp(files, 'state\d+\.mat', 'match'))/nsteps;
            end
        end
        
        %-----------------------------------------------------------------%
        function plotProgress(ensemble, range)
            % INTERNAL
            % Utility function for showing the progress of simulating a
            % range of ensemble members. Only available for
            % 'simulationType' = 'background'.
            [h_progress, h_qoi] = deal([]);
            n = 0;
            while true
                pause(0.1);
                progress = ensemble.getEnsembleMemberProgress(range);
                h_progress = plotEnsembleProgress(ensemble, progress, range, h_progress);
                if ensemble.qoi.ResultHandler.numelData > n
                    h_qoi = ensemble.plotQoI(h_qoi);
                    n = ensemble.qoi.ResultHandler.numelData;
                end
                drawnow
                if all(isinf(progress)), break; end
            end
        end
        
        %-----------------------------------------------------------------%
        function h = plotQoI(ensemble, h)
            % Creates plot(s) of the quantity of interest for all simulated
            % ensemble members
            %
            % SYNOPSIS:
            %   h = ensemble.plotQoI();
            %
            % OPTIONAL PARAMETERS:
            %   h - ??
            if nargin < 2, h = []; end
            h = ensemble.qoi.plotEnsembleQoI(ensemble, h);
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
#COPYRIGHT#
%}
