classdef TotalMobility < StateFunction
    properties
    end
    
    methods
        function gp = TotalMobility(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('Mobility');
        end
        function mobT = evaluateOnDomain(prop, model, state)
            mob = prop.getEvaluatedDependencies(state, 'Mobility');
            mobT = 0;
            for i = 1:numel(mob)
                mobT = mobT + mob{i};
            end
        end
    end
end