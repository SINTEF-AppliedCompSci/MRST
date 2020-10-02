classdef BaseSamples
    % Template class for holding the stochastic samples relevant for all
    % ensemble members in an MRSTEnsemble
    
    properties
        num = inf   % inf means that we can sample new ensemble members on the fly
        data        % precomputed sample data. Supported formats:
                    %      * Cell array of data samples
                    %      * Instance of ResultHandler class with
                    %        information about storage location and names. See
                    %        ResultHandler class and <EXAMPLE> for details
        generatorFn % Function for generating a stochastic sample
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function samples = BaseSamples(varargin)
            % Class constructor.
            samples = merge_options(samples, varargin{:});
            assert(xor(isempty(samples.data), isempty(samples.generatorFn))        , ...
                   'Please provide either a generatorFn or precomputed sample data');
            if isempty(samples.data)
                % Samples will be computed on the fly using generatorFn
            elseif iscell(samples.data)
                % We have a cell array of data samples
                samples.num = numel(samples.data);
            elseif isa(samples.data, 'ResultHandler')
                % Samples are stored to file - we are given a ResultHandler
                % for loading sample data
                samples.num = samples.data.numelData();
            else
                % Input is not on the correct format - throw an error
                error(['data must either be a cell array of samples ' , ...
                       'or an instance of the ''ResultHandler'' class']);
            end
        end
        
        %-----------------------------------------------------------------%
        function sampleData = getSample(samples, seed, problem)
            % Get a simple sample based on the seed, which is the seed for
            % the random generator or an index of data. Second input
            % argument ''problem'' is only needed if samples are computed
            % on the fly with generatorFn
            if ~isempty(samples.data)
                assert(seed >= 1 && seed <= samples.num, ...
                    ['seed must be in the range ' , ...
                     '[1, %d] = [1, samples.num]'], samples.num);
                sampleData = samples.data{seed};
            else
                sampleData = samples.generatorFn(problem, seed);
            end
        end
        
        %-----------------------------------------------------------------%
        function problem = setSample(samples, sampleData, problem) %#ok
            % Set sampleData to the problem. The specific implementation of
            % this function depends on the type of sample
            error('Template class not meant for direct use!');
        end
        
        %-----------------------------------------------------------------%
        function problem = getSampleProblem(samples, baseProblem, seed)
            % Based on the baseProblem, create a problem with sample
            % given by the seed
            sampleData = samples.getSample(seed, baseProblem);
            problem    = samples.setSample(sampleData, baseProblem);
            % Set output directory for the sample
            problem.OutputHandlers.wellSols.dataFolder = num2str(seed);
            problem.OutputHandlers.states.dataFolder   = num2str(seed);
            problem.OutputHandlers.reports.dataFolder  = num2str(seed);
            % Check if data directory exists, and make it if it does'nt
            dataDir = problem.OutputHandlers.states.getDataPath();
            if ~exist(dataDir, 'dir')
                mkdir(dataDir);
            end
        end
        
    end
end
    
