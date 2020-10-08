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
        
        simulationType = 'background';
        maxWorkers = 4;
        spmdEnsemble;
        
        verbose = false
        
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
                ensemble.directory = fullfile(mrstOutputDirectory(), 'ensemble', ensemble.setup.name);
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
        function simulateEnsembleMember(ensemble, seed, varargin)
            % Run simulation for ensemble member that corresponds to seed
            if ensemble.qoi.isComputed(seed)
                % QoI is already computed - nothing to do here!
                if ensemble.verbose
                    disp(strcat("Simulation for ", num2str(seed), " found on disc"));
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
                    ensemble.simulateEnsembleMembersBackground(range, rangePos);
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
        function simulateEnsembleMembersBackground(ensemble, range, rangePos)
%             error('Not implemented yet!');
            pathmap = capture_mrstpath();
            loadCmd = sprintf('load("%s");', fn);
            mrstPathCommand = sprintf('mrstPath(''reregister'', %s)', pathmap{:});
            for i = 1:numel(rangePos)-1
                r = range(rangePos(i):rangePos(i+1)-1);
                runCmd = sprintf('ensemble.simulateEnsemblesCore([%s]);', num2str(r));
%                 command    = sprintf('%s %s %s ', mrstPathCommand, loadCmd, runCommand);
%                 command    = sprintf('%s %s ', loadCmd, runCommand);
                runCommandsBackground({loadCmd, runCmd}, 'matlab_arg', '', 'linux_arg', '');
            end
        end
        
    end
end




