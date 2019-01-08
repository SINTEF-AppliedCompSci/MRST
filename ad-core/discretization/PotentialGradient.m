classdef PotentialGradient < GridProperty
    properties

    end
    
    methods
        function dp = evaluateOnGrid(prop, model, state)
            act = model.getActivePhases();
            nph = sum(act);
            Grad = model.operators.Grad;
            
            dp = cell(1, nph);
            if model.FlowPropertyFunctions.CapillaryPressure.pcPresent(model)
                p = model.getProp(state, 'PhasePressures');
                for i = 1:nph
                    dp{i} = Grad(p{i});
                end
            else
                p = model.getProp(state, 'pressure');
                [dp{:}] = deal(Grad(p));
            end
        end
    end
end