classdef SpatialUpscalingLevel < MLMCLevel
    
    properties
        partition
    end
    
    methods
        function level = SpatialUpscalingLevel(levelNo, partition, varargin)
            level = level@MLMCLevel(levelNo, varargin{:});
            level.constructor = @(setup) spatialUpscalingConstructor(setup, partition);
        end
    end
end