classdef WellRockSamples < WellSamples & RockSamples
    % Class that combines ensemble configurations for both rock samples and
    % well samples.
    %
    % The data property for this class is expected to be a cell array with
    % structs containing the fields 'rock' and 'well', and each of them
    % should be compatible with the data properties expected in RockSamples
    % and WellSamples, respectively.
    
    properties
        % Inherited properties only
    end
    
    methods
        
        
        %-----------------------------------------------------------------%
        function problem = setSample(samples, sampleData, problem)
            
            % Update the well settings for the problem
            problem = setSample@WellSamples(samples, sampleData.well, problem);
            
            % Update the rock of the model in the problem.
            problem = setSample@RockSamples(samples, sampleData.rock, problem);
        end
        
        
    end
end
