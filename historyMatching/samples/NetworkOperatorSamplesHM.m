classdef NetworkOperatorSamplesHM < OperatorSamplesHM & NetworkOperatorSamples
    % Class that combines ensemble configurations for both rock samples and
    % well samples.
    %
    % The data property for this class is expected to be a cell array with
    % structs containing the fields 'rock' and 'well', and each of them
    % should be compatible with the data properties expected in RockSamples
    % and WellSamples, respectively.

    
    methods
        
        
        %-----------------------------------------------------------------%
        function samples = NetworkOperatorSamplesHM(varargin)
            samples = samples@NetworkOperatorSamples(varargin{:});
        end
        
    end
end