classdef MLMCLevel

    properties
        constructor = @(setup) setup;
        solver      = @(problem) simulatePackedProblem(problem);
        levelNo
    end
    
    methods
        function level = MLMCLevel(levelNo, varargin)
            level = merge_options(level, varargin{:});
            level.levelNo = levelNo;
        end
    end
        
end