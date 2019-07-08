classdef FaceTotalMobility < StateFunction
    properties
    end
    
    methods
        function gp = FaceTotalMobility(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('FaceMobility');
        end
        function mobT = evaluateOnDomain(prop, model, state)
            mob = prop.getEvaluatedDependencies(state, 'FaceMobility');
            mobT = 0;
            for i = 1:numel(mob)
                mobT = mobT + mob{i};
            end
        end
    end
end