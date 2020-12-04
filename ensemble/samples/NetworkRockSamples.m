classdef NetworkRockSamples < RockSamples
    % Class that maps RockSamples to a network-type reduced reservoir model
    % (GPSNet, INSIM, etc). It has the same properties as RockSamples,
    % along with an additional connectionIndices property that tells us how
    % to map the values from samples.data to the different connections.
    %

    
    properties
        % Inherited properties only
        connectionIndices
    end
    
    methods
        
        
        %-----------------------------------------------------------------%
        function problem = setSample(samples, sampleData, problem)
            
            assert(~isempty(samples.connectionIndices), ...
                'Not able to map sampleData to problem since no connection indices given');
            
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