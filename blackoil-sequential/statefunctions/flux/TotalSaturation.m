classdef TotalSaturation < StateFunction
    properties

    end
    
    methods
        function gp = TotalSaturation(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('s', 'state');
        end
        function sT = evaluateOnDomain(prop, model, state)
            sT = state.sT;
        end
    end
end