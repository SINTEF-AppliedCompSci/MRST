classdef CompositeSamplesHM < BaseSamplesHM & CompositeSamples
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
    

    
    methods
        function samples = CompositeSamplesHM(parentSamples, varargin)
            samples = samples@CompositeSamples(parentSamples, varargin{:});
        end
        
        
        %-----------------------------------------------------------------%
        % Functions related to history matching
        %-----------------------------------------------------------------%
        function sampleSize = getSampleSize(samples)
            parentSizes = cellfun(@(samples) samples.getSampleSize(), samples.parentSamples);
            sampleSize = sum(parentSizes);
        end
        
        %-----------------------------------------------------------------%
        function sampleVectors = getSampleVectors(samples)
            
            sampleVectors = samples.parentSamples{1}.getSampleVectors();
            
            for i = 2:numel(samples.parentSamples)
                sampleVectors = [sampleVectors ; ...
                                 samples.parentSamples{i}.getSampleVectors()];
            end
        end
        
        %-----------------------------------------------------------------%
        function samples = setSampleVectors(samples, newSampleVectors)
            parentSizes = cellfun(@(samples) samples.getSampleSize(), samples.parentSamples);          
            sInd = [0 cumsum(parentSizes)];
            for i = 1:numel(samples.parentSamples)
                samples.parentSamples{i} = samples.parentSamples{i}.setSampleVectors(newSampleVectors((sInd(i)+1):sInd(i+1), :));
            end
            
            % Reset generator function so that the new samples are used
            fn = @(problem, seed) compositeGeneratorFn(samples, problem, seed);
            samples.generatorFn = fn;
        end
        
        
    end
    
end
