classdef PhasePotentialDifference < AutoDiffFunction
    properties

    end
    
    methods
        function gp = PhasePotentialDifference(model, varargin)
            gp@AutoDiffFunction(model, varargin{:});
            gp = gp.dependsOn('PressureGradient');
            if norm(model.gravity) > 0
                gp = gp.dependsOn('GravityPotentialDifference');
            end
        end
        function v = evaluateOnDomain(prop, model, state)
            dp = prop.getEvaluatedDependencies(state, 'PressureGradient');
            v = dp;
            if norm(model.gravity) > 0
                rhogdz = prop.getEvaluatedDependencies(state, 'GravityPotentialDifference');
                for i = 1:numel(dp)
                    v{i} = v{i} + rhogdz{i};
                end
            end
        end
    end
end