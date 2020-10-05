classdef MRSTEnsemble
    % Class that facilitates ensembles in MRST.
    %
    % SYNOPSIS:
    %   ensembles = Ensemble(baseProblem, samples, qoi)
    %
    % DESCRIPTION:
    %   Extension of `MRSTExample` class to accomodate ensembles of models.
    % 
    % REQUIRED PARAMETERS:
    %   baseProblemName - Name of function that generates the base problem,
    %                     containing all aspects that are common between
    %                     all ensemble members
    %
    %   samples         - Sample object, holding what is the specific
    %                     configurations for each individual ensemble
    %                     member, and the mapping of such a sampel onto the
    %                     baseProblem
    %
    %   qoi             - Quantity of interest object, with functionality
    %                     to read and store the relevant simulation results
    %                     from the ensemble
    %
    % 
    %
    % OPTIONAL PARAMETERS:
    %   directory       - Path to where we store results when simulating
    %                     the ensemble. Automatically generated if not
    %                     provided
    %
    %   storeOutput     - Boolean. If false, the simulation results will be
    %                     deleted after qoi is obtained.
    %
    %   solve           - Function handle for simulating/running an
    %                     ensemble member.
    %                     Default: @simulatePackedProblem(problem)
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
        qoiOptions
        
        directory % getPath() See MonteCarloSimulator.m
        storeOutput = false
        
        
        % Function handle for solving a given problem
        solve = @(problem) simulatePackedProblem(problem); 
        
        simulationType = 'background';
        maxWorkers = 4;
        spmdEnsemble;
        
    end
    
    methods
        %-----------------------------------------------------------------%
        function ensemble = MRSTEnsemble(name, samples, ...
                qoi, varargin)
            %ENSEMBLE Create an ensemble object based on a baseProblem and
            % a set of parameters/variables/properties specific to each
            % ensemble member
            
            [ensemble, extra] = merge_options(ensemble, varargin{:});
            
            if isa(name, 'MRSTExample')
                ensemble.setup = name;
            else
                ensemble.setup = MRSTExample(name, extra{:});
            end
            if isempty(ensemble.directory)
                ensemble.directory = fullfile(mrstOutputDirectory(), 'ensemble', ensemble.setup.name);
            end
%             ensemble.baseProblem = ensemble.getPackedSimulationProblem('Directory', opt.Directory, 'Name', 'baseProblem');
            
            ensemble.samples = samples;
            ensemble.num = samples.num;
            
            baseProblem = ensemble.getBaseProblem();
            ensemble.qoi = qoi.validateQoI(baseProblem);
            
            ensemble = ensemble.prepareEnsemble();
            
        end
        
        function ensemble = prepareEnsemble(ensemble)
            switch ensemble.simulationType
                case 'parallel'
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
        function setUpDataDirectory(ensemble)
            % Set up directory path, while also checking if it exists and
            % if there are already computed samples
            if isempty(ensemble.directory)
                example_hash = ensemble.getExampleHash();
                ensemble_hash = ensemble.getEnsembleHash();
                ensemble.directory = fullfile(mrstOutputDirectory(), ...
                    'ensemble', example_hash, ensemble_hash);
            end
            % Check if output directory exists
            dataPath = ensemble.getDataPath();
            if ~exist(dataPath, 'dir')
                mkdir(dataPath);
            end
        end
        
        %-----------------------------------------------------------------%
        function hash = getEnsembleHash(ensemble)
            % Get unique hash value for this Ensemble setup
            % Todo: Implement getters for hash values of stochConfigs
            % and qois that uniquely identifies their setup
            
            samplesClass = class(ensemble.samples);
            qoiClass    = class(ensemble.qoi);
            numMembers = num2str(ensemble.num);
            unhashedstr = lower(horzcat(samplesClass, qoiClass, numMembers));
            try
                % Calculate hash value
                md = java.security.MessageDigest.getInstance('SHA-256');
                hash = sprintf('%2.2x', typecast(md.digest(uint8(unhashedstr)), 'uint8')');
            catch
                % ... fallback using example string
                hash = unhashedstr;
            end
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
                    ensemble.simulateEnsembleMembersParallel(range, rangePos);
                case 'parallel'
                    ensemble.simulateEnsembleMembersParallel(range, rangePos);
                case 'background'
                    ensemble.simulateEnsembleMembersBackground(range, rangePos);
            end
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




