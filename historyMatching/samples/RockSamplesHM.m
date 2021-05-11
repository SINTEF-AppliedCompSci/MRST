classdef RockSamplesHM < BaseSamplesHM & RockSamples
    % Class for holding stochastic samples that represent uncertain rock
    % properties for an MRSTEnsemble.
    %
    % DESCRIPTION:
    %   This class is an extention of `BaseSamples` and is used within a 
    %   MRSTEnsemble to organize stochastic rock properties. These rock 
    %   properties can either be pre-computed or generated on the fly.
    % 
    % SYNOPSIS
    %   samples = RockSamples('data', data);
    %   samples = RockSamples('generatorFn', generatorFn);
    %
    % OPTIONAL PARAMETERS
    %   'data' - precomputed sample data. Supported formats:
    %               * Cell array of data samples
    %               * Instance of ResultHandler class with
    %                 information about storage location and names. See
    %                 ResultHandler class for details
    %
    %   'generatorFn' - Function for generating a stochastic sample
    % 
    % NOTE:
    %   Either 'data' or 'generatorFn' must be provided.
    %   Each sample should consist of a struct with the fields 'perm' and 
    %   'poro'.
    %
    % SEE ALSO:
    %   `WellSamples`, `DeckSamples`, `BaseSamples`, `MRSTExample`, `BaseQoI`
    
    properties
        
        minPoroValue = 0;
        minPermValue = 0;
        maxPoroValue = inf;
        maxPermValue = inf;
    end
    
    
    methods
       
        function samples = RockSamplesHM(varargin)
            samples = samples@RockSamples(varargin{:});
        end
        
        %-----------------------------------------------------------------%
        % Functions related to history matching
        %-----------------------------------------------------------------%
        function sampleVectors = getSampleVectors(samples)
            assert(~isempty(samples.data), ...
                'This function only work with pre-generated samples');
            
            exampleRock = samples.data{1};
            rockIsField = isfield(exampleRock, 'rock');
            if rockIsField
                exampleRock = samples.data{1}.rock;
            end
            
            
            assert(numel(exampleRock.poro) == numel(exampleRock.perm), ...
                'Inconsistent permeability and porosity sizes');
            
            numCells = numel(exampleRock.poro);
            sampleVectors = zeros(2*numCells, samples.num);
            
            for i = 1:samples.num
                
                if rockIsField
                    rockData = samples.data{i}.rock;
                else
                    rockData = samples.data{i};
                end
                
                if samples.transformSampleVectors
                    % Logarithmic transform of permeability 
                    perm = convertTo(rockData.perm(:), milli*darcy);
                    sampleVectors(1:numCells, i) = log(perm);
                else
                    sampleVectors(1:numCells, i) = rockData.perm(:);
                end
                sampleVectors(numCells+1:numCells*2, i) = rockData.poro(:);
            end
        end
        
        %-----------------------------------------------------------------%
        function samples = setSampleVectors(samples, newSampleVectors)
            
            assert(size(newSampleVectors, 2) == samples.num, ...
                'number of columns of new samples does not match ensemble size');
            
            rockIsField = isfield(samples.data{1}, 'rock');
            
            if rockIsField
                numCells = numel(samples.data{1}.rock.poro);
            else
                numCells = numel(samples.data{1}.poro);
            end
            
            assert(size(newSampleVectors, 1) == numCells*2, ...
                'number of rows of new sample does not match old sample size');
            
            for i = 1:samples.num
                
                if samples.transformSampleVectors
                    perm = exp(newSampleVectors(1:numCells, i));
                    permData = convertFrom(perm, milli*darcy);
                else
                    permData = newSampleVectors(1:numCells, i);
                end
                                
                poroData = newSampleVectors(numCells+1:numCells*2, i);

                % Check and fix potentially bad values
                permData(permData < samples.minPermValue) = samples.minPermValue;
                permData(permData > samples.maxPermValue) = samples.maxPermValue;
                poroData(poroData < samples.minPoroValue) = samples.minPoroValue;
                poroData(poroData > samples.maxPoroValue) = samples.maxPoroValue;
                
                if rockIsField
                    samples.data{i}.rock.perm(:) = permData;
                    samples.data{i}.rock.poro(:) = poroData;
                else
                    samples.data{i}.perm(:) = permData;
                    samples.data{i}.poro(:) = poroData;
                end
            end
            
        end

        
    end
    
end

%{
#COPYRIGHT#
%}