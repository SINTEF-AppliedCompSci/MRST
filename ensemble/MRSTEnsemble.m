classdef MRSTEnsemble < MRSTExample
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
        baseProblem
        num
        samples
        qoi
        qoiOptions
        
        directory % getPath() See MonteCarloSimulator.m
        storeOutput = false
        
        
        % Function handle for solving a given problem
        solve = @(problem) simulatePackedProblem(problem); 
        
        
        parallel = true;
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
            
            opt = struct('Directory', []);
            [opt, extra] = merge_options(opt, varargin{:});
            
            ensemble = ensemble@MRSTExample(name, extra{:});
            if isempty(opt.Directory)
                opt.Directory = fullfile(mrstOutputDirectory(), 'ensemble', ensemble.name);
            end
            ensemble.baseProblem = ensemble.getPackedSimulationProblem('Directory', opt.Directory, 'Name', 'baseProblem');
            
            ensemble.samples = samples;
            ensemble.num = samples.num;
            
            ensemble.qoi = qoi.validateQoI(ensemble.baseProblem);
            
            if ensemble.parallel
                spmd
                    spmdEnsemble = ensemble;
                end
                ensemble.spmdEnsemble = spmdEnsemble;
            end
            
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
            dataPath = ensemble.baseProblem.OutputHandlers.states.dataDirectory();
        end
               
        %-----------------------------------------------------------------%
        function simulateEnsembleMember(ensemble, seed, varargin)
            % Run simulation for ensemble member that corresponds to seed
            if ensemble.qoi.isComputed(seed)
                % QoI is already computed - nothing to do here!
                return;
            end
            % Set up sample problem from seed
            problem = ensemble.samples.getSampleProblem(ensemble.baseProblem, seed);
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
            if ensemble.parallel
                ensemble.simulateEnsembleMembersParallel(range, rangePos);
            else
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
        function simulateEnsembleMembersParallel(ensemble, range, rangePos)
            spmdEns = ensemble.spmdEnsemble;
            spmd
                spmdRange = range(rangePos(labindex):rangePos(labindex+1)-1);
                spmdEns.simulateEnsembleMembersCore(spmdRange);
            end
        end
        
        %-----------------------------------------------------------------%
        function simulateEnsembleMembersBackground(ensemble, range, rangePos)
            error('Not implemented yet!');
%             fn = fullfile(ensemble.getDataPath(), 'ensemble.mat');
%             if ~exist(fn, 'file')
%                 save(fn, 'ensemble');
%             end
%             pathmap = capture_mrstpath();
%             loadCommand = sprintf('load("%s");', fn);
%             mrstPathCommand = sprintf('mrstPath(''reregister'', %s)', pathmap{:});
%             for i = 1:numel(rangePos)-1
%                 r = range(rangePos(i):rangePos(i+1)-1);
%                 runCommand = sprintf('ensemble.simulateEnsemblesCore([%s]);', num2str(r));
%                 command    = sprintf('%s %s %s ', mrstPathCommand, loadCommand, runCommand);
%                 command    = sprintf('%s %s ', loadCommand, runCommand);
%                 runCommandsBackground(command, 'matlab_arg', '', 'linux_arg', '');
%             end
        end
        
    end
end

function pathmap = capture_mrstpath()
    mods  = mrstPath();
    paths = mrstPath(mods{:});

    pathmap = [ reshape(mods , 1, []) ; ...
                reshape(paths, 1, []) ];
end


