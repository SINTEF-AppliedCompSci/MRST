classdef PotentialGradient < GridProperty
    properties

    end
    
    methods
        function v = evaluateOnGrid(prop, model, state)
            gdz = model.getGravityGradient();
            rho = model.getProp(state, 'Density');
        end
    end
end