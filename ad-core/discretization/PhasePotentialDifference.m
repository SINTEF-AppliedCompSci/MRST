classdef PhasePotentialDifference < GridProperty
    properties

    end
    
    methods
        function v = evaluateOnGrid(prop, model, state)
            dp = model.getProp(state, 'PressureGradient');
            v = dp;
            if norm(model.gravity) > 0
                rhogdz = model.getProps(state, 'GravityPotentialDifference');
                for i = 1:numel(dp)
                    v{i} = v{i} + rhogdz{i};
                end
            end
        end
    end
end