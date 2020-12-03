classdef CompositeSamples < BaseSamples
    % Class for composite samples from a number of parent samples.
    %
    %
    % DESCRIPTION:
    %   This class is used within a MRSTEnsemble to organize composite
    %   stochastic samples that are made up of multiple samples, e.g., rock
    %   samples and oil-water contact samples
    %
    % EXAMPLE:
    % See `ensembleExampleDoubleSampleSets`
    %
    % SEE ALSO:
    %   `RockSamples`, `WellSamples`, `DeckSamples`, `MRSTExample`, `BaseQoI`
    
    properties
        parentSamples         % Cell array of samples making up the 
                              % composite sample
        tensorProduct = false % If false, the same seed is used to generate
                              % all parent samples. This requires that the
                              % parentSamples.num are equal. If true, use
                              % all combinations of all parent samples.
                              % Seed will be interpreted as a linear index
                              % into the multiindex space spanned by all
                              % possible tensor products of
                              % (1:n1) x (1:n2) x ... x (1:nm),
                              % where ni =parentSamples{i}.num
                              % It is up to the user to ensure that this 
                              % makes sense, e.g., that samples are 
                              % stochastically independent.
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function samples = CompositeSamples(parentSamples, varargin)
            samples.parentSamples = parentSamples;
            samples = merge_options(samples, varargin{:});
            % Get number of samples
            num = samples.numParentSamples;
            if ~all(cellfun(@(samples) isempty(samples.generatorFn), samples.parentSamples))
                    samples.tensorProduct = false;
            end
            if samples.tensorProduct
                samples.num = prod(num);
            else
                % Seed will be interpreted as the combined sample from each
                % parent generated with the same seed
                assert(all(num == num(1)), ...
                    ['For non-tensor composite samples, all '            , ...
                     'parentSamples must have the same number of samples'])
                samples.num = num(1);
            end
            % Set generator function
            fn = @(problem, seed) compositeGeneratorFn(samples, problem, seed);
            samples.generatorFn = fn;
        end
        
        %-----------------------------------------------------------------%
        function problem = setSample(samples, sampleData, problem)
            % Set sample by looping through all parents
            for i = 1:samples.numParents
                data    = sampleData.([class(samples.parentSamples{i}), 'Data']);
                problem = samples.parentSamples{i}.setSample(data, problem);
            end
        end
        
        %-----------------------------------------------------------------%
        % Utilities for getting parentSample information
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function n = numParents(samples)
            n = numel(samples.parentSamples);
        end
        
        %-----------------------------------------------------------------%
        function n = numParentSamples(samples)
            n = cellfun(@(samples) samples.num, samples.parentSamples);
        end
        
        
      
    end
    
    
    
    methods(Access=protected)
        %-------------------------------------------------------------------------%
        function sampleData = compositeGeneratorFn(samples, problem, seed)
            % Function for getting a composite sample
            sampleData = struct();
            if samples.tensorProduct
                % Create tensor product mapping and find parentSeeds
                n      = num2cell(samples.numParentSamples);
                n      = cellfun(@(n) 1:n, n, 'UniformOutput', false);
                [n{:}] = ndgrid(n{:}); 
                n      = cellfun(@(n) n(:), n, 'UniformOutput', false);
                parentSeed = cellfun(@(n) n(seed), n);
            else
                % Parent seeds are equal to seed
                parentSeed = repmat(seed, samples.numParents, 1);
            end
            for i = 1:samples.numParents
                % Store parent sample data with field name [class, 'Data']
                parentData = samples.parentSamples{i}.getSample(parentSeed(i), problem);
                sampleData.([class(samples.parentSamples{i}), 'Data']) = parentData;
            end
        end
    end
end

