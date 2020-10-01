classdef MRSTEnsemble < MRSTExample
    % Class that facilitates ensembles in MRST.
    %
    % SYNOPSIS:
    %   ensembles = Ensemble(baseProblem, configurations, qoi)
    %
    % DESCRIPTION:
    %   Extension of `PhysicalModel` class to accomodate reservoir-specific
    %   features such as fluid and rock as well as commonly used phases and
    %   variables.
    %
    % REQUIRED PARAMETERS:
    %   baseProblemName - Name of function that generates the base problem,
    %                     containing all aspects that are common between
    %                     all ensemble members
    %
    %   configurations  - Configuration object, holding the the specific
    %                     configurations for the individual ensemble
    %                     parameters, and the mapping of such a
    %                     configuration to the baseProblem
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
        stochasticConfigs
        qoi
        qoiOptions
        
        directory % getPath() See MonteCarloSimulator.m
        storeOutput = false
        
        
        % Function handle for solving a given problem
        solve = @(problem) simulatePackedProblem(problem); 
        
        
        
        
    end
    
    methods
        %-----------------------------------------------------------------%
        function ensemble = MRSTEnsemble(baseProblemName, stochasticConfigs, ...
                qoi, varargin)
            %ENSEMBLE Create an ensemble object based on a baseProblem and
            % a set of parameters/variables/properties specific to each
            % ensemble member
            
            ensemble = ensemble@MRSTExample(baseProblemName, varargin{:});

            ensemble.baseProblem = ensemble.getPackedSimulationProblem();
            
            ensemble.stochasticConfigs = stochasticConfigs;
            ensemble.num = stochasticConfigs.num;
            
            ensemble.qoi = qoi.validateQoI(ensemble.baseProblem);
            
            ensemble.setUpDataDirectory();
            
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
            
            stochConfigClass = class(ensemble.stochasticConfigs);
            qoiClass    = class(ensemble.qoi);
            numMembers = num2str(ensemble.num);
            unhashedstr = lower(horzcat(stochConfigClass, qoiClass, numMembers));
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
            dataPath = ensemble.directory;
        end
        
        %-----------------------------------------------------------------%
        function ok = saveQOI(ensemble, problem, seed, qoi)
            error('not implemented');
            
            % construct the directory for the given problem's qoi
            qoi_path = 'something';
            ensemble.qoi.save(qoi)
            
        end
               
        %-----------------------------------------------------------------%
        function simulateEnsembleMember(ensemble, seed, vargin)
            % run simulation for ensemble member that corresponds to seed
            
            % TODO: Check if QoI exists
            %     if so - return
            
            
            opt = struct('clearPackedSimulationOutput', 'true');
            [opt, extra] = merge_options(opt, vargin{:});
            
            problem = ensemble.stochasticConfigs.getProblem(...
                ensemble.baseProblem, seed);
            ensemble.solve(problem);
            
            currentQoI = ensemble.qoi.getQOI(problem);
            ensemble.qoi.save(currentQoI);
            
            if opt.clearPackedSimulationOutput
                clearPackedSimulatorOutput(problem, 'prompt', false);
            end
            
        end
    end
end



