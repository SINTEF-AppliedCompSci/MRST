classdef MRSTEnsemble
    % Class that facilitates ensembles in MRST.
    %
    % SYNOPSIS:
    %   ensembles = Ensemble(baseProblem, samples, qoi)
    %
    % DESCRIPTION:
    %   This class is used to organize, run and post-process ensemble
    %   simulations in MRST. It can be used for a varity of problems, such
    %   as uncertainty quantification or history matching.
    % 
    % REQUIRED PARAMETERS:
    %   mrstExample - Instance of MRSTExample, or a name that generates a
    %                 MRSTExample, that serves as a baseProblem for the
    %                 ensemble, containing all aspects that are common
    %                 between all ensemble members.
    %
    %   samples     - Sample object, holding what is the specific
    %                 configurations for each individual ensemble
    %                 member, and the mapping of such a sampel onto the
    %                 baseProblem
    %
    %   qoi         - Quantity of interest object, with functionality
    %                 to read and store the relevant simulation results
    %                 from the ensemble
    %
    % 
    %
    % OPTIONAL PARAMETERS:
    %   directory        - Path to where we store results when simulating
    %                      the ensemble. Automatically generated if not
    %                      provided
    %                      Default: mrstOutputDirectory()/ensemble/baseProblemName
    %
    %   storeOutput      - Boolean. If false, the simulation results will be
    %                      deleted after qoi is obtained.
    %                      Default: false
    %
    %   solve            - Function handle for simulating/running an
    %                      ensemble member.
    %                      Default: @simulatePackedProblem(problem)
    %
    %   simulationType   - String to control how to run the ensemble.
    %                      Acceptable values:
    %                        'serial'     - No parallelization
    %                        'parallel'   - Run ensemble in parallel using 
    %                            the Parallel Computing Toolbox
    %                        'background' - Run ensemble members by spawning
    %                            new matlab sessions in the background
    %
    %   maxWorkers       - Maximum number of parallel workers to use for
    %                      processing the ensemble members
    %
    %   verbose          - Boolean (default: false). Print some extra info
    %
    %   deleteOldResults - Boolean (default: false). Delete any old
    %                      simulation results that might exist.
    %
    %   extra parameters that are passed on to the MRSTExample class, in
    %   case mrstExample is a string.
    %
    % RETURNS:
    %   Class instance.
    %
    % SEE ALSO:
    %   `MRSTExample`, 
    properties
        setup
        num
        samples
        qoi
        
        directory % getPath() See MonteCarloSimulator.m
        storeOutput = false
        deleteOldResults = false
        
        
        % Function handle for solving a given problem
        solve = @(problem) simulatePackedProblem(problem); 
        
        simulationType = 'serial';
        maxWorkers = 4;
        spmdEnsemble;
        
        verbose = false
        
        evalFn       = @simulateEnsembleMembersStandalone
        matlabBinary = ''
        plotProgress
        
    end
    
    methods
        %-----------------------------------------------------------------%
        function ensemble = MRSTEnsemble(mrstExample, samples, qoi, varargin)
            %ENSEMBLE Create an ensemble object based on a baseProblem and
            % a set of parameters/variables/properties specific to each
            % ensemble member
            
            [ensemble, extra] = merge_options(ensemble, varargin{:});
            
            if isa(mrstExample, 'MRSTExample')
                ensemble.setup = mrstExample;
            else
                ensemble.setup = MRSTExample(mrstExample, extra{:});
            end
            
            % Set up directory
            if isempty(ensemble.directory)
                ensemble.directory = ensemble.getDefaultPath();
            end
            if ensemble.deleteOldResults && exist(ensemble.directory, 'dir')
                rmdir(ensemble.directory, 's');
            end
            
            
            ensemble.samples = samples;
            ensemble.num = samples.num;
            
            baseProblem = ensemble.getBaseProblem();
            ensemble.qoi = qoi.validateQoI(baseProblem);
            
            ensemble = ensemble.prepareEnsemble();
            
        end
        
        function ensemble = prepareEnsemble(ensemble)
            switch ensemble.simulationType
                case 'parallel'
                    % Check if we have started a parallel session already,
                    % and whether it has the correct number of workers
                    if isempty(gcp('nocreate'))
                        parpool(ensemble.maxWorkers);
                    elseif gcp('nocreate').NumWorkers ~= ensemble.maxWorkers
                        delete(gcp);
                        parpool(ensemble.maxWorkers);
                    end
                    
                    spmd
                        spmdEns = ensemble;
                    end
                    ensemble.spmdEnsemble = spmdEns;
                case 'background'
                    fn = fullfile(ensemble.getDataPath(), 'ensemble.mat');
                    if ~exist(fn, 'file')
                        save(fn, 'ensemble');
                    end
            end
        end
        
        function problem = getBaseProblem(ensemble)
            problem = ensemble.setup.getPackedSimulationProblem('Directory', ensemble.directory, 'Name', 'baseProblem');
        end
        
        %-----------------------------------------------------------------%
        function memberProblem = getProblem(ensemble, i)
            % Get the specific problem with configurations for ensemble 
            % member i
            assert(i > 0 && i <= ensemble.num, 'illegal ensemble ID');
            memberProblem = ensemble.configurations.getProblem(i);
        end
        
        %-----------------------------------------------------------------%
        function ensemble = setSolver(ensemble, solve)
            ensemble.solve = solve;
        end
        
        %-----------------------------------------------------------------%
        function dataPath = getDataPath(ensemble)
            dataPath = ensemble.directory();
        end
              
        %-----------------------------------------------------------------%
        function defaultPath = getDefaultPath(ensemble)
            defaultPath = fullfile(mrstOutputDirectory(), 'ensemble', ensemble.setup.name);
        end
        
        %-----------------------------------------------------------------%
        function simulateEnsembleMember(ensemble, seed, varargin)
            % Run simulation for ensemble member that corresponds to seed
            if ensemble.qoi.isComputed(seed)
                % QoI is already computed - nothing to do here!
                if ensemble.verbose
                    disp(strcat("Simulation for ", num2str(seed), " found on disk"));
                end
                return;
            end
            % Get base problem
            baseProblem = ensemble.getBaseProblem();
            % Set up sample problem from seed
            problem = ensemble.samples.getSampleProblem(baseProblem, seed);
            % Solve problem
            ensemble.solve(problem);
            % Compute QoI
            ensemble.qoi.getQoI(problem);
            
            if ~ensemble.storeOutput
                clearPackedSimulatorOutput(problem, 'prompt', false);
            end
            
        end
        
        %-----------------------------------------------------------------%
        function simulateEnsembleMembers(ensemble, range, varargin)
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
                    ensemble.simulateEnsembleMembersSerial(range, rangePos);
                case 'parallel'
                    ensemble.simulateEnsembleMembersParallel(range, rangePos);
                case 'background'
                    ensemble.simulateEnsembleMembersBackground(range, rangePos, varargin{:});
            end
        end
        
        %-----------------------------------------------------------------%
        function simulateAllEnsembleMembers(ensemble, varargin)
            assert(~isinf(ensemble.num), ...
                'Ensemble size not define, please use ensemble.simulateEnsembleMembers(range) instead');
            ensemble.simulateEnsembleMembers(1:ensemble.num, varargin);
        end
        
        %-----------------------------------------------------------------%
        function simulateEnsembleMembersCore(ensemble, range)
           for seed = reshape(range, 1, [])
               ensemble.simulateEnsembleMember(seed);
           end
        end

        %-----------------------------------------------------------------%
        function simulateEnsembleMembersSerial(ensemble, range, rangePos)
            ensemble.simulateEnsembleMembersCore(range);
        end
        
        %-----------------------------------------------------------------%
        function simulateEnsembleMembersParallel(ensemble, range, rangePos)
            spmdEns = ensemble.spmdEnsemble;
            spmd
                spmdRange = range(rangePos(labindex):rangePos(labindex+1)-1);
                spmdEns.simulateEnsembleMembersCore(spmdRange);
            end
        end
        
        %-----------------------------------------------------------------%
        function simulateEnsembleMembersBackground(ensemble, range, rangePos, varargin)
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
                evalFunWrapper(ensemble.evalFn, {fileName, r}         , ...
                               'progressFileNm', progressFileNm(r)    , ...
                               'moduleList'    , moduleList           , ...
                               'matlabBinary'  , ensemble.matlabBinary);
            end
            fprintf(['Started %d new Matlab sessions. ', ...
                     'Waiting for simulations ...\n'  ], n);
            pause(0.1);
            if opt.plotProgress
                nr = numel(num2str(max(range)));
                em = ['%', num2str(nr), 'u'];    
                while true
                    progress = ensemble.getEnsembleMemberProgress(range);
                    % TODO: Implement function for plotting progress 
                    % plotEnsembleSimulationProgress(h, ensemble, progress)
                    clc;
                    fprintf('Simulating ensemble members ... \n')
                    fprintf(repmat('-', 1, numel('Ensemble member ') + nr + 3 + 100 + 10))
                    fprintf('\n');
                    for i = 1:numel(range)
                        p = min(progress(i), 1);
                        fprintf(['Ensemble member ', em, ': |'], range(i));
                        fprintf(repmat('=', 1, round(p*100))    );
                        fprintf(repmat(' ', 1, round((1-p)*100)));
                        fprintf('| %6.2f %% \n', p*100);
                    end
                    fprintf(repmat('-', 1, 100 + numel('Ensemble member ') + nr + 3 + 10))
                    fprintf('\n');
                    pause(0.5)
                    if all(isinf(progress)), break; end
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function progress = getEnsembleMemberProgress(ensemble, range)
            if nargin < 2, range = ensemble.num; end
            progress = zeros(numel(range),1);
            nsteps   = numel(ensemble.setup.schedule.step.val);
            for i = 1:numel(range)
                if exist(fullfile(ensemble.directory(), ...
                                  ['qoi', num2str(range(i)), '.mat']), 'file')
                    progress(i) = inf;
                    continue
                end
                dataDir = fullfile(ensemble.directory(), num2str(range(i)));
                if ~exist(dataDir, 'dir')
                    continue;
                end
                files = ls(dataDir);
                progress(i) = numel(regexp(files, 'state\d+\.mat', 'match'))/nsteps;
            end
        end
        
    end
end




