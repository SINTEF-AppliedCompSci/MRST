classdef TotalSaturation < StateFunction
    properties

    end
    
    methods
        function gp = TotalSaturation(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('s', 'state');
        end
        function sT = evaluateOnDomain(prop, model, state)
            sT = 0;
            s = model.getProp(state, 's');
            for i = 1:numel(s)
                sT = sT + s{i};
            end
        end
    end
end