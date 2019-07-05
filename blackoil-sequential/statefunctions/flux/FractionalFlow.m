classdef FractionalFlow < StateFunction
    properties
    end
    
    methods
        function gp = FractionalFlow(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'Mobility', 'TotalMobility'}, 'FlowPropertyFunctions');
        end
        function f = evaluateOnDomain(prop, model, state)
            [mob, mobT] = prop.getEvaluatedDependencies(state, 'Mobility');
            f = mob;
            for i = 1:numel(mob)
                f{i} = f{i}./mobT;
            end
        end
    end
end