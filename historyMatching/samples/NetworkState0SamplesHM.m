classdef NetworkState0SamplesHM < State0SamplesHM & NetworkState0Samples
    % Class that combines ensemble configurations for both rock samples and
    % well samples.
    %
    % The data property for this class is expected to be a cell array with
    % structs containing the fields 'rock' and 'well', and each of them
    % should be compatible with the data properties expected in RockSamples
    % and WellSamples, respectively.

    
    methods
        
        
        %-----------------------------------------------------------------%
        function samples = NetworkState0SamplesHM(varargin)
            samples = samples@NetworkState0Samples(varargin{:});
        end
        
    end
end