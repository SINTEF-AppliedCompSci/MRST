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
        
        minPoroValue = eps;
        minPermValue = eps;
        maxPoroValue = inf;
        maxPermValue = inf;
        
        useKozenyCarman = false;
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
            
            if samples.useKozenyCarman
                sampleVectors = sampleVectors(numCells+1:numCells*2, :);
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
            
            if samples.useKozenyCarman
                assert(size(newSampleVectors, 1) == numCells, ...
                    'number of rows of new sample does not match old sample size');
                
                % Expand vector sizes
                newSampleVectors = repmat(newSampleVectors, 2, 1);
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
                
                if samples.useKozenyCarman
                    permData =  poroData.^3.*(1e-5)^2./(0.81*72*(1-poroData).^2);
                end
                
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
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
