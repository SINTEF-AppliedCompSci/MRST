classdef FaceTotalMobility < StateFunction
    properties (Access = private)
        mobility_name = 'FaceMobility';
    end
    
    methods
        function gp = FaceTotalMobility(model, mobility_name)
            if nargin < 2
                mobility_name = 'FaceMobility';
            end
            gp@StateFunction(model);
            gp = gp.dependsOn(mobility_name);
            gp.mobility_name = mobility_name;
        end
        function mobT = evaluateOnDomain(prop, model, state)
            mob = prop.getEvaluatedDependencies(state, prop.mobility_name);
            mobT = 0;
            for i = 1:numel(mob)
                mobT = mobT + mob{i};
            end
        end
    end
end