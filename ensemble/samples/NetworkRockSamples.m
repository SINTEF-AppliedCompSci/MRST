classdef NetworkRockSamples < RockSamples
    % Class that combines ensemble configurations for both rock samples and
    % well samples.
    %
    % The data property for this class is expected to be a cell array with
    % structs containing the fields 'rock' and 'well', and each of them
    % should be compatible with the data properties expected in RockSamples
    % and WellSamples, respectively.
    
    properties
        % Inherited properties only
        connectionIndices
    end
    
    methods
        
        
        %-----------------------------------------------------------------%
        function problem = setSample(samples, sampleData, problem)
            
            numConnections = numel(samples.connectionIndices.cells);
            assert(numel(sampleData.poro) == numConnections, ...
                'mismatch between number of connections and sample data size');
            
            sampleRock = problem.SimulatorSetup.model.rock;
            
            for c = 1:numConnections
                sampleRock.poro(samples.connectionIndices.cells{c}) = sampleData.poro(c);
                sampleRock.perm(samples.connectionIndices.cells{c}) = sampleData.perm(c);
            end
            
            % Update the rock of the model in the problem.
            problem = setSample@RockSamples(samples, sampleRock, problem);
        end
        
        
    end
end